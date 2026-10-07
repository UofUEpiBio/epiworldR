#ifndef EPIWORLD_RUN_OUTPUTS_HPP
#define EPIWORLD_RUN_OUTPUTS_HPP

/**
 * @brief A column-oriented table, the in-memory form of one output CSV.
 *
 * Columns have the same names, order and rows as the matching CSV file
 * (see `run_output_names()`), so a binding can turn a table into a data
 * frame, a dict of arrays, etc. with a single converter.
 */
struct OutputTable
{
    using Column = std::variant<
        std::vector< int >,
        std::vector< double >,
        std::vector< std::string >
    >;

    std::vector< std::string > colnames;
    std::vector< Column > columns;

    template<typename T>
    void add(std::string name, std::vector< T > values)
    {
        colnames.push_back(std::move(name));
        columns.emplace_back(std::move(values));
    }

    size_t nrow() const
    {
        return columns.empty() ?
            0u :
            std::visit([](const auto & c) { return c.size(); }, columns[0u]);
    }
};

/// Tables of one simulation (or, in `SaverMemory::results()`, of all of
/// them), keyed by output name.
using RunOutputs = std::map< std::string, OutputTable >;

/// The names of the outputs, which are also the suffixes of the CSV files.
inline const std::vector< std::string > & run_output_names()
{
    static const std::vector< std::string > names = {
        "virus_info", "virus_hist", "tool_info", "tool_hist", "total_hist",
        "transmission", "transition", "reproductive", "generation",
        "active_cases", "outbreak_size", "hospitalizations"
    };
    return names;
}

/**
 * @brief Writes a table as a space-separated file with a header.
 *
 * Strings are double-quoted, except for the `*_sequence` columns.
 *
 * In `EPI_DEBUG` builds, each line is prefixed with the id of the thread that
 * wrote it (a leading `thread` column), as the CSVs always were. It describes
 * the writer, not the simulation, so it is not part of `OutputTable`.
 */
inline void write_table(const std::string & fn, const OutputTable & table)
{

    std::ofstream file(fn, std::ios_base::out);
    if (!file)
        throw std::runtime_error("Could not open file \"" + fn + "\" for writing.");

    #ifdef EPI_DEBUG
    file << "thread ";
    #endif

    for (size_t j = 0u; j < table.colnames.size(); ++j)
        file << table.colnames[j] << (j + 1u < table.colnames.size() ? " " : "\n");

    // Quoting only the string columns that are not sequences
    std::vector< bool > quote(table.columns.size());
    for (size_t j = 0u; j < quote.size(); ++j)
        quote[j] =
            std::holds_alternative< std::vector< std::string > >(table.columns[j]) &&
            table.colnames[j].find("_sequence") == std::string::npos;

    for (size_t i = 0u; i < table.nrow(); ++i)
    {

        #ifdef EPI_DEBUG
        file << EPI_GET_THREAD_ID() << " ";
        #endif

        for (size_t j = 0u; j < table.columns.size(); ++j)
        {
            std::visit(
                [&](const auto & col) { file << (quote[j] ? "\"" : "") << col[i] << (quote[j] ? "\"" : ""); },
                table.columns[j]
            );
            file << (j + 1u < table.columns.size() ? " " : "\n");
        }

    }

}

/**
 * @brief Keeps the outputs of `run_multiple()` in memory.
 *
 * An object of this class is a valid `fun` for `Model::run_multiple()`. Each
 * simulation's tables are stored by `sim_id`, so results do not depend on
 * which thread finishes first. Copies share the same storage (which is what
 * lets `std::function` and OpenMP copy it); `run_multiple()` already
 * serializes calls to `fun`, but the object is not otherwise thread-safe.
 *
 * @code{.cpp}
 * SaverMemory saver({"total_hist", "transition"});
 * model.run_multiple(100, 50, 123, saver, true, true, 4);
 * RunOutputs all = saver.results();
 * @endcode
 */
class SaverMemory
{
private:

    struct State
    {
        std::vector< std::string > whats;
        std::vector< RunOutputs > runs; ///< Indexed by sim_id
    };

    std::shared_ptr< State > state = std::make_shared< State >();

public:

    /// @param whats Names of the outputs to keep (see `run_output_names()`).
    SaverMemory(std::vector< std::string > whats = {"total_hist"})
    {
        state->whats = std::move(whats);
    }

    template<typename TSeq>
    void operator()(size_t sim_id, Model<TSeq> * model) const
    {
        auto out = model->get_db().get_run_outputs(state->whats);
        if (sim_id >= state->runs.size())
            state->runs.resize(sim_id + 1u);
        state->runs[sim_id] = std::move(out);
    }

    /// All simulations, concatenated, with a leading 0-based `sim_id` column.
    RunOutputs results() const
    {

        RunOutputs ans;

        for (const auto & w : state->whats)
        {

            OutputTable * all = nullptr;

            for (size_t s = 0u; s < state->runs.size(); ++s)
            {

                auto it = state->runs[s].find(w);
                if (it == state->runs[s].end())
                    continue;

                const OutputTable & t = it->second;

                if (all == nullptr)
                {
                    all = &(ans[w] = t);
                    all->colnames.insert(all->colnames.begin(), "sim_id");
                    all->columns.insert(
                        all->columns.begin(),
                        OutputTable::Column(std::vector< int >(t.nrow(), static_cast<int>(s)))
                    );
                    continue;
                }

                std::get< std::vector< int > >(all->columns[0u]).insert(
                    std::get< std::vector< int > >(all->columns[0u]).end(),
                    t.nrow(), static_cast<int>(s)
                );

                for (size_t j = 0u; j < t.columns.size(); ++j)
                    std::visit(
                        [&](auto & dst)
                        {
                            const auto & src = std::get< std::decay_t< decltype(dst) > >(t.columns[j]);
                            dst.insert(dst.end(), src.begin(), src.end());
                        },
                        all->columns[j + 1u]
                    );

            }

        }

        return ans;

    }

    /// Forgets the stored simulations (call it before reusing the saver).
    void clear() { state->runs.clear(); }

};

#endif
