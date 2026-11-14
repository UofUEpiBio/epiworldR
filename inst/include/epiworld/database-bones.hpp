#ifndef EPIWORLD_DATABASE_BONES_HPP
#define EPIWORLD_DATABASE_BONES_HPP

class Model;

/// Every disease: the default of the database's `disease` arguments (see
/// `Model::disease_id()`).
constexpr int all_diseases = -1;

/**
 * @brief Which coordinates a transition matrix is over (see
 * `DataBase::get_transition_probability()`).
 * @details `pair`: status × stage, as the engine records it. `stage` and
 * `status`: the pair's counts summed over the other coordinate, so a change
 * of status within a stage counts as staying in the stage matrix (and a
 * change of stage within a status, in the status matrix).
 */
enum class Margin { pair, stage, status };

/**
 * @brief A transition matrix over a margin, with its labels.
 * @details Row `from`, column `to`: `values[from + to * size()]`. The labels
 * of row and column `i` are `status[i]` and `stage[i]`; `stage` is empty in a
 * status matrix and `status` in a stage matrix.
 */
struct TransitionMatrix {
    Margin margin = Margin::pair;
    std::string disease;             ///< The disease ("" for `Margin::status`)
    bool normalized = true;          ///< Probabilities (rows sum to 1) or counts
    std::vector< std::string > status;
    std::vector< std::string > stage;
    std::vector< epiworld_double > values;

    size_t size() const { return margin == Margin::status ? status.size() : stage.size(); }
    epiworld_double operator()(size_t from, size_t to) const { return values.at(from + to * size()); }
    std::string label(size_t i) const; ///< "status / stage", the stage, or the status
    void print() const;
};

class Virus;

class UserData;

/**
 * @brief The variants and the tools, by id.
 * @details Kept once per model: a variant or a tool is not tied to one
 * disease's records.
 */
struct DataRegistry {
    // Variants, by id: the viruses added to the model first (founders, kept
    // even if removed, so ids stay put), then the run's mutations
    size_t n_founders = 0u;
    std::vector< std::string > virus_name;
    std::vector< int > virus_origin_date;
    std::vector< int > virus_parent_id;

    // Tools, by id
    std::vector< std::string > tool_name;
    std::vector< int > tool_origin_date;

    bool operator==(const DataRegistry & other) const;
};

/**
 * @brief What is recorded about one disease: counts by state (today and by
 * day), the variants' and tools' histories, transitions, transmissions and
 * hospitalizations.
 * @details One per model for now; with several diseases (NH-5b), one per
 * disease.
 */
struct DiseaseRecord {
    HospitalizationsTracker m_hospitalizations; ///< Tracker for hospitalization events

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_virus;

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_tool;

    // {Susceptible, Infected, etc.}
    std::vector< int > today_total;

    // Totals
    int today_total_nviruses_active = 0;

    // Variants history
    std::vector< int > hist_virus_date;
    std::vector< int > hist_virus_id;
    std::vector< epiworld_fast_uint > hist_virus_state;
    std::vector< int > hist_virus_counts;

    // Tools history
    std::vector< int > hist_tool_date;
    std::vector< int > hist_tool_id;
    std::vector< epiworld_fast_uint > hist_tool_state;
    std::vector< int > hist_tool_counts;

    // Overall hist
    std::vector< int > hist_total_date;
    std::vector< int > hist_total_nviruses_active;
    std::vector< epiworld_fast_uint > hist_total_state;
    std::vector< int > hist_total_counts;
    std::vector< int > hist_transition_matrix;

    // Transmission network
    std::vector< int > transmission_date;                 ///< Date of the transmission event
    std::vector< int > transmission_source;               ///< Id of the source
    std::vector< int > transmission_target;               ///< Id of the target
    std::vector< int > transmission_virus;              ///< Id of the variant
    std::vector< int > transmission_source_exposure_date; ///< Date when the source acquired the variant
    std::vector< int > transmission_layer;                ///< Id of the layer of the contact (-1: none, e.g., seeded)

    std::vector< int > transition_matrix;

    /**
     * @name Declared transitions (natural history)
     * @details How many agents took each declared transition, since the last
     * record (see `get_hist_transitions()`).
     */
    ///@{
    std::vector< int > transition_today;       ///< [transition]
    std::vector< int > hist_transition_date;
    std::vector< int > hist_transition_id;
    std::vector< int > hist_transition_counts;
    ///@}

    bool operator==(const DiseaseRecord & other) const;
};

/**
 * @brief Statistical data about the process
 *
 * **Implementation details:**
 * <a href="../impl/data-collection-and-analysis.md">Data Collection and Analysis</a>,
 * <a href="../impl/generation-interval-and-reproductive-number.md">Generation Interval and Reproductive Number Calculation</a>
 *
 * @details The database is the one entry point (`Model::get_db()`). It
 * holds its data in parts: the variants and tools (`DataRegistry`) and the
 * disease's records (`DiseaseRecord`).
 */
class DataBase {
    friend class Model;
private:
    Model * model;

    DataRegistry reg;  ///< Variants and tools
    DiseaseRecord rec; ///< The disease's records

    int sampling_freq = 1;

    UserData user_data;

    void update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo = false
    );

    void update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo = false
    );

    void update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo = false
    );

    void record_transition(epiworld_fast_uint from, epiworld_fast_uint to, bool undo);

    /// One more agent took the declared transition `id` (see
    /// `get_hist_transitions()`)
    void record_declared_transition(int id) { ++rec.transition_today[static_cast< size_t >(id)]; }

    /**
     * @brief The status and the stage of a state of the engine (the one place
     * the pair is split; models without statuses have the status "Free").
     */
    void split_state(
        epiworld_fast_uint state,
        std::string & status,
        std::string & stage
    ) const;

    /// `split_state()` over a vector of states.
    void split_states(
        const std::vector< epiworld_fast_uint > & states,
        std::vector< std::string > * status,
        std::vector< std::string > * stage
    ) const;

    /// Today's total of a stage (`by_stage`) or of a status.
    int today_total_of(const std::string & name, bool by_stage) const;

    /// Checks a disease argument: `all_diseases`, or a disease of the model.
    void check_disease(int disease) const;

    /// The transmissions, written into arrays (see `get_transmissions()`).
    void get_transmissions_into(
        int * date,
        int * source,
        int * target,
        int * virus,
        int * source_exposure_date,
        int * layer
    ) const;

    /// The disease column of `n` rows (if `names` is not `nullptr`).
    void disease_rows(std::vector< std::string > * names, size_t n) const;


public:

    #ifdef EPI_DEBUG
    int n_transmissions_potential = 0;
    int n_transmissions_today     = 0;
    #endif

    DataBase() = delete;
    DataBase(Model & m) : model(&m), user_data(m) {};
    DataBase(const DataBase & db);
    // DataBase & operator=(const DataBase & m);

    /**
     * @brief Records a new variant, giving the virus its id
     *
     * @param v The virus.
     * @param founder `true` when `Model::add_virus()` records the virus: it
     * starts its own lineage (and the variants of the last run are dropped).
     * Otherwise `v` has mutated: the new variant is a child of the id it had,
     * and its agent's counts move from the parent to it.
     */
    void record_virus(Virus & v, bool founder = false);
    void record_tool(Tool & t); ///< Gives the tool the next id.
    void reset();
    void drop_variants(); ///< Keeps the founders' records, dropping the run's variants.
    Model * get_model();
    void record();

    const std::vector< int > & get_nexposed() const;
    size_t size() const;

    /**
     * @name Get recorded information from the model
     *
     * @details The counts are by status and stage: one row per pair (and
     * day, and virus or tool, in the histories). Models without statuses
     * report the status "Free". A total by stage or by status is a sum over
     * the rows.
     *
     * @param stage A stage's name.
     * @param status A status's name.
     * @return In `get_today_total_stage(stage)` and
     * `get_today_total_status(status)`, how many agents are in that stage (in
     * any status) or in that status (in any stage) now. A stage and a status
     * may share a name; each function looks only among its own.
     * @return In `get_today_total`, the current counts of each pair.
     * @return In `get_today_virus`, the current counts of each pair for
     * each virus.
     * @return In `get_hist_total`, the time series of each pair.
     * @return In `get_hist_virus` and `get_hist_tool`, the time series of
     * each pair for each virus and tool.
     * @return In `get_active_cases`, the number of active cases (currently infected individuals)
     * for each virus at each point in time.
     *
     * @param disease Which disease's rows (`Model::disease_id()`), or
     * `all_diseases` (the default) for every disease's. An index that is not
     * a disease of the model is an error.
     * @param disease_names If not `nullptr`, the disease of each row (its
     * name), one per row of the other outputs. The same two arguments end
     * the other functions that return rows (transitions, transmissions,
     * generation times, hospitalizations).
     */
    ///@{
    int get_today_total_stage(
        const std::string & stage,
        int disease = all_diseases
    ) const;
    int get_today_total_status(const std::string & status) const;
    void get_today_total(
        std::vector< std::string > * status,
        std::vector< std::string > * stage,
        std::vector< int > * counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_today_virus(
        std::vector< std::string > & status,
        std::vector< std::string > & stage,
        std::vector< int > & id,
        std::vector< int > & counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_hist_total(
        std::vector< int > * date,
        std::vector< std::string > * status,
        std::vector< std::string > * stage,
        std::vector< int > * counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_hist_virus(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & status,
        std::vector< std::string > & stage,
        std::vector< int > & counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_hist_tool(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & status,
        std::vector< std::string > & stage,
        std::vector< int > & counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    /**
     * @brief How many agents took each declared transition of the natural
     * history, by day.
     * @details One row per recorded day and transition, the transitions in
     * the order they were declared. Seeding and initial placements are not
     * transitions. Requires a natural history.
     */
    void get_hist_transitions(
        std::vector< int > & date,
        std::vector< std::string > & transition,
        std::vector< int > & counts,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_hist_transition_matrix(
        std::vector< std::string > & status_from,
        std::vector< std::string > & stage_from,
        std::vector< std::string > & status_to,
        std::vector< std::string > & stage_to,
        std::vector< int > & date,
        std::vector< int > & counts,
        bool skip_zeros,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_active_cases(
        std::vector< int > & date,
        std::vector< int > & virus_id,
        std::vector< int > & count,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    void get_outbreak_size(
        std::vector< int > & date,
        std::vector< int > & virus_id,
        std::vector< int > & size,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;
    ///@}

    /**
     * @brief Get the transmissions object
     *
     * @param date
     * @param source
     * @param target
     * @param virus
     * @param source_exposure_date
     * @param layer Id of the layer of the contact that transmitted the virus
     * (-1 for infections that did not come from a contact, e.g., seeds)
     */
    ///@{
    void get_transmissions(
        std::vector<int> & date,
        std::vector<int> & source,
        std::vector<int> & target,
        std::vector<int> & virus,
        std::vector<int> & source_exposure_date,
        std::vector<int> & layer,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

    ///@}

    void write_data(
        std::string fn_virus_info,
        std::string fn_virus_hist,
        std::string fn_tool_info,
        std::string fn_tool_hist,
        std::string fn_total_hist,
        std::string fn_transmission,
        std::string fn_transition,
        std::string fn_reproductive_number,
        std::string fn_generation_time,
        std::string fn_active_cases,
        std::string fn_outbreak_size,
        std::string fn_hospitalizations
        ) const;

    /**
     * @brief In-memory tables of this simulation, one per requested output
     * @param whats Names of the outputs (see `run_output_names()`).
     * @details Each table has exactly the columns and rows of the CSV that
     * `write_data()` writes for it. Throws `std::invalid_argument` on an
     * unknown name.
     */
    RunOutputs get_run_outputs(const std::vector< std::string > & whats) const;

    /***
     * @brief Record a transmission event
     * @param i,j Integers. Id of the source and target agents.
     * @param virus Integer. Id of the virus.
     * @param i_expo_date Integer. Date when the source agent was infected.
     * @details
     * If i is -1, then it means that the agent was assigned a virus at the
     * beginning of the simulation.
     */
    void record_transmission(int i, int j, int virus, int i_expo_date, int layer = -1);

    size_t get_n_viruses() const; ///< Get the number of viruses
    size_t get_n_tools() const; ///< Get the number of tools

    void set_user_data(std::vector< std::string > names);
    void add_user_data(std::vector< epiworld_double > x);
    void add_user_data(epiworld_fast_uint j, epiworld_double x);
    UserData & get_user_data();


    /**
     * @brief Computes the reproductive number of each case
     *
     * @details By definition, whereas it computes R0 (basic reproductive number)
     * or Rt/R (the effective reproductive number) will depend on whether the
     * virus is allowed to circulate naïvely or not, respectively.
     *
     * @param fn File where to write out the reproductive number.
     * @details
     * In the case of `MapVec_type<int,int>`, the key is a vector of 3 integers:
     * - Virus id
     * - Source id
     * - Date when the source was infected
     * 
     * The seed cases are those with source id != -1 and date when the source
     * was infected == 0. Cases with source id == -1 is equivalent to the 
     * Model's reproductive number, this is, the initial number of cases
     * at the beginning of the simulation.
     */
    ///@{
    MapVec_type<int,int> get_reproductive_number(int disease = all_diseases) const;

    void get_reproductive_number(
        std::string fn
        ) const;
    ///@}

    /**
     * @brief The transitions over the whole run, as probabilities or counts
     * @param margin `Margin::pair` (status × stage), `Margin::stage` or
     * `Margin::status`.
     * @param disease The disease (`Model::disease_id()`). `all_diseases`, the
     * default, is that one disease when the model has one; with several, a
     * pair or stage matrix needs the disease (a status matrix does not: it is
     * the same for every disease).
     * @param normalize Probabilities (each row with transitions sums to 1),
     * or the counts.
     * @details Entry (i, j) counts the agent-days that started in i and ended
     * in j (staying counts, on the diagonal).
     */
    TransitionMatrix get_transition_probability(
        Margin margin = Margin::pair,
        int disease = all_diseases,
        bool normalize = true
    ) const;

    bool operator==(const DataBase & other) const;
    bool operator!=(const DataBase & other) const {return !operator==(other);};

    /**
     * Calculates the generating time
     * @param agent_id,virus_id,time,gentime vectors where to save the values
     *
     * @details
     * The generation time is the time between the infection of the source and
     * the infection of the target.
    */
   ///@{
    void get_generation_time(
        std::vector< int > & agent_id,
        std::vector< int > & virus_id,
        std::vector< int > & time,
        std::vector< int > & gentime,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const; ///< Get the generation time

    void get_generation_time(
        std::string fn
    ) const; ///< Write the generation time to a file
    ///@}

    /**
     * @brief Record a hospitalization event for an agent.
     * 
     * @param agent Reference to the agent being hospitalized.
     * 
     * @details
     * For each hospitalization, the method records:
     * - The current date from the model
     * - The virus ID from the agent's virus
     * - For each tool the agent has, a separate record with weight = 1/N
     *   where N is the number of tools
     * - If the agent has no tools, a single record with tool_id = -1 and
     *   weight = 1.0
     */
    void record_hospitalization(Agent & agent);

    /**
     * @brief Get the full time series of hospitalization data.
     * 
     * @param date Output vector for dates.
     * @param virus_id Output vector for virus IDs.
     * @param tool_id Output vector for tool IDs.
     * @param count Output vector for counts (number of hospitalized individuals).
     * @param weight Output vector for summed weights (fractional contribution
     *   based on tool distribution).
     * 
     * @details
     * Returns the full time series of hospitalization data. For each unique 
     * (virus_id, tool_id) combination observed, returns an entry for every 
     * day from 0 to ndays-1.
     * 
     * The `count` vector contains the actual number of individuals hospitalized 
     * for that (date, virus_id, tool_id) combination. This is useful for 
     * answering questions like "how many total people were hospitalized?" 
     * regardless of their tools.
     * 
     * The `weight` vector contains fractional contributions: if an agent has N 
     * tools, each tool gets weight = 1/N. Summing weights across all tool_ids 
     * for a given date and virus_id gives the total number of hospitalizations.
     */
    void get_hospitalizations(
        std::vector<int> & date,
        std::vector<int> & virus_id,
        std::vector<int> & tool_id,
        std::vector<int> & count,
        std::vector<double> & weight,
        int disease = all_diseases,
        std::vector< std::string > * disease_names = nullptr
    ) const;

};


#endif
