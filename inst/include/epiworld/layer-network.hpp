#ifndef EPIWORLD_LAYER_NETWORK_HPP
#define EPIWORLD_LAYER_NETWORK_HPP

namespace layer {

/**
 * @brief A fixed network layer (see `Network`).
 * @param name Name of the layer (unique in a model).
 * @param graph The ties; its number of vertices must match the model's
 * population. Directed if the graph is.
 * @ingroup contact_layers
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
std::shared_ptr< Network<TSeq> > network(std::string name, const AdjList & graph);

/**
 * @brief A network layer drawn at every reset (see `Network`).
 * @param name Name of the layer (unique in a model).
 * @param generator Called with the model at every reset; the graph it returns
 * must have one vertex per agent.
 * @ingroup contact_layers
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
std::shared_ptr< Network<TSeq> > network(
    std::string name,
    typename Network<TSeq>::Generator generator
);

/**
 * @brief Contact layer with explicit ties (a network).
 *
 * @details The contacts of an agent in this layer are its ties. A network
 * layer is created with `layer::network()`, either from a graph, which is then
 * *fixed*, or from a *generator*, a function that receives the model and
 * returns a graph, which the layer calls at every reset: each run (and each
 * replicate of `run_multiple()`) then gets its own draw, from its own seed.
 * The `graph::` functions return generators; any function
 * `AdjList(Model<TSeq> &)` is one.
 *
 * ### Ties and their order
 *
 * Each agent's ties are kept in the order they were added (a graph is loaded
 * with each agent's ties in ascending id order, without repeats), and that
 * order is load-bearing: the single draw that picks the infector walks the
 * contacts in this order, so whichever agent sits at a given position is the
 * one recorded. Adding or removing a tie never permutes the others. Agents with
 * more than `EPI_NEIGHBOR_INDEX_THRESHOLD` ties also get a position index, so
 * that looking a tie up is not a scan.
 *
 * ### Edits
 *
 * `add_edge()`, `rm_edge()` and the rewiring function are the only way to
 * change the ties, and they keep the model's queue and agents-by-state index in
 * step, so they are safe at any point of a run. An edit made during a run
 * applies to that run only: the next reset restores a fixed graph, and redraws
 * a generated one. An edit made between runs changes a fixed graph for the runs
 * that follow; a generated graph is redrawn anyway.
 *
 * ### Rewiring
 *
 * `set_rewiring()` gives the layer a function that is called at the end of
 * every step (not on day 0), after the global events, with the layer and the
 * model; `graph::rewire_degseq()` returns the degree-preserving one.
 *
 * @tparam TSeq Type of sequence.
 * @ingroup contact_layers
 */
template<typename TSeq>
class Network : public ContactLayer<TSeq>
{
    friend class Model<TSeq>;
    friend class Queue<TSeq>;

public:

    /// Returns the graph of a run (called at every reset).
    using Generator = std::function< AdjList(Model<TSeq> &) >;

    /// Changes the ties at the end of each step.
    using Rewiring = std::function< void(Network<TSeq> &, Model<TSeq> &) >;

private:

    Generator generator; ///< Empty for a fixed graph
    Rewiring rewiring;
    bool directed = false;

    std::vector< std::vector< size_t > > ties; ///< Neighbors of each agent

    /// Position of each neighbor, only for agents with many ties
    std::unordered_map< size_t, std::unordered_map< size_t, size_t > > index;

    // The fixed graph the next reset restores, one row per agent
    std::vector< size_t > baseline_start;
    std::vector< size_t > baseline_ids;
    bool baseline_stale = false; ///< Edited between runs

    void set_graph(const AdjList & graph);
    void capture_baseline();
    void restore_baseline();
    void build_index(size_t i);
    size_t find(size_t i, size_t j) const; ///< Position of `j` in `ties[i]`
    bool append(size_t i, size_t j, bool check);
    void erase_at(size_t i, size_t pos);
    void check_endpoints(size_t i, size_t j) const;
    void edit();

protected:

    // Created with layer::network()
    Network(std::string name, const AdjList & graph);
    Network(std::string name, Generator generator);

    template<typename T>
    friend std::shared_ptr< Network<T> > network(std::string, const AdjList &);

    template<typename T>
    friend std::shared_ptr< Network<T> > network(
        std::string, typename Network<T>::Generator
    );

public:

    size_t size() const { return ties.size(); } ///< Number of agents
    bool is_directed() const { return directed; }

    /// Number of ties of agent `i` (in a directed network, the ties it keeps).
    size_t get_degree(size_t i) const { return ties[i].size(); }

    /**
     * @brief The neighbors of agent `i`, without allocating.
     * @details In order (see the class description). The view is invalidated by
     * any change to the agent's ties. The layer must belong to a model.
     */
    NeighborsView<TSeq> neighbors_view(size_t i);

    /**
     * @name Edit the ties
     *
     * @details Undirected ties only, changed at both ends: `add_edge()` and
     * `rm_edge()` throw in a directed network.
     *
     * @return Whether a tie was added (removed).
     * @throws std::range_error if an id is out of range.
     * @throws std::logic_error if `i == j` or the network is directed.
     */
    ///@{
    bool add_edge(size_t i, size_t j);
    bool rm_edge(size_t i, size_t j);
    ///@}

    /// Whether `i` and `j` are tied (in a directed network, whether `i -> j`).
    bool has_edge(size_t i, size_t j) const;

    /**
     * @brief Swaps two ties: `j`, `i`'s neighbor at position `pos_i`, and `l`,
     * `k`'s neighbor at position `pos_k`, trade places.
     *
     * @details (i-j), (k-l) become (i-l), (k-j); in an undirected network the
     * other ends follow. No degree changes. Used by rewiring functions.
     * @throws std::range_error if an id or a position is out of range.
     * @throws std::logic_error if the network is undirected and the ties share
     * an endpoint (`i`, `j`, `k` and `l` must be four different agents).
     */
    void swap_ties(size_t i, size_t pos_i, size_t k, size_t pos_k);

    /// Sets the function that rewires the ties at the end of each step.
    void set_rewiring(Rewiring fun) { rewiring = std::move(fun); }
    bool has_rewiring() const { return static_cast< bool >(rewiring); }

    /**
     * @brief The ties as an edge list.
     * @details Each undirected tie once (`source <= target`); each directed
     * tie as kept by its source. The vectors are appended to.
     */
    void get_edgelist(
        std::vector< int > & source,
        std::vector< int > & target
    ) const;

    /// Writes the edge list (see `get_edgelist()`) to a file.
    void write_edgelist(const std::string & fn) const;

    /**
     * @brief The neighbors of `agent` that carry a virus, written to `out`
     * (room for `get_degree(agent)` entries).
     * @details Neighbors in a state flagged in `exclude` (if not null) are
     * skipped. With `Record`, each contact is reported to the model's
     * post-sampling callback.
     * @return The number of contacts written.
     */
    template< bool Record >
    size_t pull(
        size_t agent,
        Model<TSeq> & model,
        const std::vector< bool > * exclude,
        size_t * out
    ) const;

    /// The network's part of a pull (see `ContactLayer::sample()`).
    size_t sample(
        size_t agent,
        Model<TSeq> & model,
        const std::vector< bool > * exclude,
        std::vector< size_t > & out
    ) override
    {

        const size_t first = out.size();
        out.resize(first + ties[agent].size());

        // The recording path is chosen once per call, not per contact
        const size_t n = model.post_sampling_on ?
            pull<true>(agent, model, exclude, out.data() + first) :
            pull<false>(agent, model, exclude, out.data() + first);

        out.resize(first + n);
        return n;

    }

    void reset(Model<TSeq> & model) override;
    void update(Model<TSeq> & model) override;

    layer::Network<TSeq> * as_network() override { return this; }

    /// An undirected network: a carrier's ties are its contacts' ties too.
    bool supports_queue() const override { return !directed; }

    std::shared_ptr< ContactLayer<TSeq> > clone() const override
    {
        return std::make_shared< Network<TSeq> >(*this);
    }

};

template<typename TSeq>
inline std::shared_ptr< Network<TSeq> > network(
    std::string name,
    const AdjList & graph
)
{
    return std::shared_ptr< Network<TSeq> >(
        new Network<TSeq>(std::move(name), graph)
    );
}

template<typename TSeq>
inline std::shared_ptr< Network<TSeq> > network(
    std::string name,
    typename Network<TSeq>::Generator generator
)
{
    return std::shared_ptr< Network<TSeq> >(
        new Network<TSeq>(std::move(name), std::move(generator))
    );
}

template<typename TSeq>
inline Network<TSeq>::Network(std::string name, const AdjList & graph) :
    ContactLayer<TSeq>(std::move(name))
{
    set_graph(graph);
    capture_baseline();
}

template<typename TSeq>
inline Network<TSeq>::Network(std::string name, Generator generator_) :
    ContactLayer<TSeq>(std::move(name)),
    generator(std::move(generator_))
{
    if (!generator)
        throw std::invalid_argument(
            "Network layer '" + this->get_name() + "': empty generator."
        );
}

template<typename TSeq>
inline void Network<TSeq>::set_graph(const AdjList & graph)
{

    const size_t n = graph.vcount();

    // AdjList::is_directed() throws on a graph with no vertices
    directed = (n > 0u) && graph.is_directed();

    const auto & dat = graph.get_dat();

    // Counting sort of the tie ends by agent: row i (from start[i] to
    // start[i + 1]) holds i's neighbors. An undirected tie i - j is written to
    // both rows, so the result does not depend on the list being symmetric.
    std::vector< size_t > start(n + 1u, 0u);
    for (size_t i = 0u; i < n; ++i)
    {
        start[i + 1u] += dat[i].size();
        if (!directed)
            for (const auto & link : dat[i])
                start[static_cast< size_t >(link.first) + 1u]++;
    }

    std::partial_sum(start.begin(), start.end(), start.begin());

    std::vector< size_t > ids(start[n]);
    std::vector< size_t > next(start.begin(), start.end() - 1);
    for (size_t i = 0u; i < n; ++i)
        for (const auto & link : dat[i])
        {
            size_t j = static_cast< size_t >(link.first);
            ids[next[i]++] = j;
            if (!directed)
                ids[next[j]++] = i;
        }

    // Each row becomes its agent's ties: sorted and without repeats
    ties.assign(n, {});
    index.clear();
    for (size_t i = 0u; i < n; ++i)
    {

        auto first = ids.begin() + static_cast< std::ptrdiff_t >(start[i]);
        auto last  = ids.begin() + static_cast< std::ptrdiff_t >(start[i + 1u]);

        std::sort(first, last);
        last = std::unique(first, last);

        ties[i].assign(first, last);

        if (ties[i].size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
            build_index(i);

    }

}

template<typename TSeq>
inline void Network<TSeq>::capture_baseline()
{

    const size_t n = ties.size();
    baseline_start.assign(n + 1u, 0u);
    for (size_t i = 0u; i < n; ++i)
        baseline_start[i + 1u] = baseline_start[i] + ties[i].size();

    baseline_ids.resize(baseline_start[n]);
    for (size_t i = 0u; i < n; ++i)
        std::copy(
            ties[i].begin(), ties[i].end(),
            baseline_ids.begin() + static_cast< std::ptrdiff_t >(baseline_start[i])
        );

}

template<typename TSeq>
inline void Network<TSeq>::restore_baseline()
{

    const size_t n = baseline_start.size() - 1u;
    ties.resize(n);
    index.clear();
    for (size_t i = 0u; i < n; ++i)
    {

        ties[i].assign(
            baseline_ids.begin() + static_cast< std::ptrdiff_t >(baseline_start[i]),
            baseline_ids.begin() + static_cast< std::ptrdiff_t >(baseline_start[i + 1u])
        );

        if (ties[i].size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
            build_index(i);

    }

}

template<typename TSeq>
inline void Network<TSeq>::build_index(size_t i)
{

    auto & pos = index[i];
    pos.clear();
    pos.reserve(ties[i].size());
    for (size_t k = 0u; k < ties[i].size(); ++k)
        pos[ties[i][k]] = k;

}

template<typename TSeq>
inline size_t Network<TSeq>::find(size_t i, size_t j) const
{

    const auto & row = ties[i];

    if (row.size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
    {
        const auto & pos = index.at(i);
        auto it = pos.find(j);
        return (it == pos.end()) ? row.size() : it->second;
    }

    for (size_t k = 0u; k < row.size(); ++k)
        if (row[k] == j)
            return k;

    return row.size();

}

template<typename TSeq>
inline bool Network<TSeq>::append(size_t i, size_t j, bool check)
{

    if (check && (find(i, j) != ties[i].size()))
        return false;

    ties[i].push_back(j);

    if (ties[i].size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
    {
        auto it = index.find(i);
        if (it == index.end())
            build_index(i);
        else
            it->second[j] = ties[i].size() - 1u;
    }

    return true;

}

template<typename TSeq>
inline void Network<TSeq>::erase_at(size_t i, size_t pos)
{

    auto & row = ties[i];
    row.erase(row.begin() + static_cast< std::ptrdiff_t >(pos));

    // Everything after `pos` shifted down by one, so the index has to follow.
    // Moving the last tie into the hole would be O(1) but would permute the
    // others, and their order decides which infector the draw picks.
    if (row.size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
        build_index(i);
    else
        index.erase(i);

}

template<typename TSeq>
inline void Network<TSeq>::check_endpoints(size_t i, size_t j) const
{

    if ((i >= ties.size()) || (j >= ties.size()))
        throw std::range_error(
            "Network layer '" + this->get_name() + "': agent ids must be below " +
            std::to_string(ties.size()) + "; got " + std::to_string(i) +
            " and " + std::to_string(j) + "."
        );

    if (i == j)
        throw std::logic_error(
            "Network layer '" + this->get_name() + "': an agent cannot be tied "
            "to itself (agent " + std::to_string(i) + ")."
        );

    if (directed)
        throw std::logic_error(
            "Network layer '" + this->get_name() + "': add_edge() and rm_edge() "
            "change both ends of a tie, which is not meaningful in a directed "
            "network."
        );

}

template<typename TSeq>
inline void Network<TSeq>::edit()
{

    // Between runs, a fixed graph is edited from its baseline: the edits of the
    // last run are dropped first, and the next reset takes the result as the
    // new baseline.
    const bool dirty = this->is_dirty();
    if (!this->edits_baseline() || generator)
        return;

    if (dirty)
    {
        restore_baseline();
        this->mark_clean();
    }

    baseline_stale = true;

}

template<typename TSeq>
inline NeighborsView<TSeq> Network<TSeq>::neighbors_view(size_t i)
{

    if (this->model == nullptr)
        throw std::logic_error(
            "Network layer '" + this->get_name() + "': neighbors_view() needs "
            "the layer to be added to a model."
        );

    if (i >= ties.size())
        throw std::range_error(
            "Network layer '" + this->get_name() + "': agent " +
            std::to_string(i) + " is out of range."
        );

    const auto & row = ties[i];
    return NeighborsView<TSeq>(
        row.empty() ? nullptr : row.data(), row.size(),
        &this->model->get_agents()
    );

}

template<typename TSeq>
inline bool Network<TSeq>::add_edge(size_t i, size_t j)
{

    check_endpoints(i, j);
    edit();

    // Two statements, so that the second end is visited even when the first
    // already had the tie
    const bool here  = append(i, j, true);
    const bool there = append(j, i, true);

    if (!here && !there)
        return false;

    if (this->model != nullptr)
        this->model->network_tie_changed(i, here ? 1 : 0, j, there ? 1 : 0);

    return true;

}

template<typename TSeq>
inline bool Network<TSeq>::rm_edge(size_t i, size_t j)
{

    check_endpoints(i, j);
    edit();

    const size_t pos_i = find(i, j);
    const size_t pos_j = find(j, i);
    const bool rm_i = pos_i != ties[i].size();
    const bool rm_j = pos_j != ties[j].size();

    // Nothing to unwind if the two were never tied
    if (!rm_i && !rm_j)
        return false;

    if (rm_i)
        erase_at(i, pos_i);

    if (rm_j)
        erase_at(j, pos_j);

    if (this->model != nullptr)
        this->model->network_tie_changed(i, rm_i ? -1 : 0, j, rm_j ? -1 : 0);

    return true;

}

template<typename TSeq>
inline bool Network<TSeq>::has_edge(size_t i, size_t j) const
{

    if ((i >= ties.size()) || (j >= ties.size()))
        throw std::range_error(
            "Network layer '" + this->get_name() + "': agent ids must be below " +
            std::to_string(ties.size()) + "; got " + std::to_string(i) +
            " and " + std::to_string(j) + "."
        );

    return find(i, j) != ties[i].size();

}

template<typename TSeq>
inline void Network<TSeq>::swap_ties(
    size_t i,
    size_t pos_i,
    size_t k,
    size_t pos_k
)
{

    if ((i >= ties.size()) || (k >= ties.size()))
        throw std::range_error(
            "Network layer '" + this->get_name() + "': agent ids must be below " +
            std::to_string(ties.size()) + "."
        );

    if (pos_i >= ties[i].size())
        throw std::range_error(
            "The tie you want to swap is out of range. Agent " +
            std::to_string(i) + " only has " + std::to_string(ties[i].size()) +
            " ties."
        );

    if (pos_k >= ties[k].size())
        throw std::range_error(
            "The tie you want to swap is out of range. Agent " +
            std::to_string(k) + " only has " + std::to_string(ties[k].size()) +
            " ties."
        );

    const size_t j = ties[i][pos_i];
    const size_t l = ties[k][pos_k];

    // The other ends follow in an undirected network, so ties that share an
    // agent would write into the rows being swapped (and leave self-ties)
    if (!directed && ((i == k) || (j == l) || (i == l) || (k == j)))
        throw std::logic_error(
            "Network layer '" + this->get_name() + "': swap_ties() needs four "
            "different agents in an undirected network; the ties (" +
            std::to_string(i) + "-" + std::to_string(j) + ") and (" +
            std::to_string(k) + "-" + std::to_string(l) + ") share one."
        );

    edit();

    std::swap(ties[i][pos_i], ties[k][pos_k]);

    if (!directed)
    {

        // The two agents that were swapped away have to be told: the one that
        // pointed back at `i` now points at `k`, and vice versa. Positions are
        // looked up rather than cached, so ties can be edited without a
        // back-pointer table going stale.
        const size_t back_j = find(j, i);
        const size_t back_l = find(l, k);

        #ifdef EPI_DEBUG
        if ((back_j >= ties[j].size()) || (back_l >= ties[l].size()))
            throw std::logic_error(
                "[epi-debug] swap_ties: the tie is not reciprocated.");
        #endif

        ties[j][back_j] = k;
        ties[l][back_l] = i;

        // Four lists may have changed. Rebuilding is O(degree), the same order
        // as the swap itself, and sidesteps the aliasing cases (an agent on
        // both sides of the swap) that piecemeal updates get wrong.
        for (size_t a : {i, k, j, l})
            if (ties[a].size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
                build_index(a);

    }
    else
    {

        for (size_t a : {i, k})
            if (ties[a].size() > EPI_NEIGHBOR_INDEX_THRESHOLD)
                build_index(a);

    }

    // Nobody's degree changed; only the queue needs to follow
    if (this->model != nullptr)
        this->model->network_ties_swapped(i, j, k, l, directed);

}

template<typename TSeq>
inline void Network<TSeq>::get_edgelist(
    std::vector< int > & source,
    std::vector< int > & target
) const
{

    for (size_t i = 0u; i < ties.size(); ++i)
        for (size_t j : ties[i])
            if (directed || (i <= j))
            {
                source.push_back(static_cast< int >(i));
                target.push_back(static_cast< int >(j));
            }

}

template<typename TSeq>
inline void Network<TSeq>::write_edgelist(const std::string & fn) const
{

    std::ofstream efile(fn, std::ios_base::out);
    efile << "source target\n";

    for (size_t i = 0u; i < ties.size(); ++i)
        for (size_t j : ties[i])
            if (directed || (i <= j))
                efile << i << " " << j << "\n";

}

template<typename TSeq>
template< bool Record >
inline size_t Network<TSeq>::pull(
    size_t agent,
    Model<TSeq> & model,
    const std::vector< bool > * exclude,
    size_t * out
) const
{

    size_t n = 0u;
    const auto & row = ties[agent];

    // During a run, most neighbors do not carry a virus. The model's compact
    // per-agent flag rules them out reading one byte each instead of loading
    // the neighbor's `Agent`; the contacts kept, and their order, are the same.
    const char * carrier = model.get_agent_carrier();
    auto & agents = model.get_agents();

    for (size_t j : row)
    {

        if ((carrier != nullptr) && (carrier[j] == 0))
            continue;

        auto & neighbor = agents[j];

        if ((exclude != nullptr) && (*exclude)[neighbor.get_state()])
            continue;

        if (neighbor.get_virus() == nullptr)
            continue;

        // The post-sampling callback sees every eligible contact
        if constexpr (Record)
            model.register_sampled_contact(j, agent);

        out[n++] = j;

    }

    return n;

}

template<typename TSeq>
inline void Network<TSeq>::reset(Model<TSeq> & model)
{

    if (generator)
        set_graph(generator(model));
    else if (this->is_dirty())
        restore_baseline();
    else if (baseline_stale)
        capture_baseline(); // Edits made between runs become the baseline

    baseline_stale = false;

    if (ties.size() != model.size())
        throw std::length_error(
            "Network layer '" + this->get_name() + "': the graph has " +
            std::to_string(ties.size()) + " vertices, but the model has " +
            std::to_string(model.size()) + " agents."
        );

    this->mark_clean();

}

template<typename TSeq>
inline void Network<TSeq>::update(Model<TSeq> & model)
{

    // Day 0 is recorded, not stepped: rewiring starts at the end of day 1
    if (rewiring && (model.today() > 0))
        rewiring(*this, model);

}

}

/**
 * @brief Rewires a share of a network layer's ties, keeping every degree.
 *
 * @details Picks two agents with probability proportional to their degree and
 * one tie of each, uniformly, and swaps them (`layer::Network::swap_ties()`),
 * skipping swaps that would create self-ties or repeated ties;
 * `floor(proportion * sum of degrees)` attempts. Use it through
 * `graph::rewire_degseq()`.
 */
template<typename TSeq>
inline void rewire_degseq(
    layer::Network<TSeq> & net,
    Model<TSeq> & model,
    epiworld_double proportion
)
{

    const size_t n = net.size();

    #ifdef EPI_DEBUG
    std::vector< size_t > degree0(n);
    for (size_t i = 0u; i < n; ++i)
        degree0[i] = net.get_degree(i);
    #endif

    // Agents with ties, weighted by degree
    std::vector< size_t > non_isolates;
    std::vector< epiworld_double > weights;
    epiworld_double nedges = 0.0;

    for (size_t i = 0u; i < n; ++i)
    {
        if (net.get_degree(i) > 0u)
        {
            non_isolates.push_back(i);
            epiworld_double wtemp = static_cast< epiworld_double >(net.get_degree(i));
            weights.push_back(wtemp);
            nedges += wtemp;
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (size_t i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    const size_t N = non_isolates.size();
    int nrewires = floor(proportion * nedges);
    while (nrewires-- > 0)
    {

        // Picking egos
        epiworld_double prob = model.runif();
        int id0 = static_cast< int >(N) - 1;
        for (size_t i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id0 = static_cast< int >(i);
                break;
            }

        prob = model.runif();
        int id1 = static_cast< int >(N) - 1;
        for (size_t i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id1 = static_cast< int >(i);
                break;
            }

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= static_cast< int >(N))
            id1 = 0;

        const size_t a0 = non_isolates[id0];
        const size_t a1 = non_isolates[id1];

        // Picking alters (relative location in their lists), uniformly
        const size_t id01 = model.runif_index(net.get_degree(a0));
        const size_t id11 = model.runif_index(net.get_degree(a1));

        const size_t alter0 = net.neighbors_view(a0).id(id01);
        const size_t alter1 = net.neighbors_view(a1).id(id11);

        // After the swap a0 is tied to alter1 and a1 to alter0. Skip swaps of
        // the same alter, and those that would tie an agent to itself or
        // repeat a tie.
        if ((alter0 == alter1) || (alter0 == a1) || (alter1 == a0))
            continue;

        if (net.has_edge(a0, alter1) || net.has_edge(a1, alter0))
            continue;

        net.swap_ties(a0, id01, a1, id11);

    }

    #ifdef EPI_DEBUG
    for (size_t i = 0u; i < n; ++i)
        if (degree0[i] != net.get_degree(i))
            throw std::logic_error("[epi-debug] Degree does not match after rewire_degseq.");
    #endif

}

namespace graph {

/**
 * @name Network generators
 *
 * @details Each returns a `layer::Network<TSeq>::Generator`: a function that
 * draws the graph for a model, with one vertex per agent, using the model's
 * random numbers. The arguments are structure, used once per draw (see the
 * contact-layers plan, point 17).
 */
///@{

/// Small-world network: ring lattice with `k` ties per agent, rewired with
/// probability `p` (see `rgraph_smallworld()`).
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline typename layer::Network<TSeq>::Generator smallworld(
    epiworld_fast_uint k,
    epiworld_double p,
    bool directed = false
)
{
    return [k, p, directed](Model<TSeq> & model) -> AdjList {
        return rgraph_smallworld(
            static_cast< epiworld_fast_uint >(model.size()), k, p, directed, model
        );
    };
}

/// Bernoulli (Erdős–Rényi) network: each tie with probability `p` (see
/// `rgraph_bernoulli()`).
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline typename layer::Network<TSeq>::Generator bernoulli(
    epiworld_double p,
    bool directed = false
)
{
    return [p, directed](Model<TSeq> & model) -> AdjList {
        return rgraph_bernoulli(
            static_cast< epiworld_fast_uint >(model.size()), p, directed, model
        );
    };
}

/// Stochastic block model (see `rgraph_sbm()`); the block sizes must add up to
/// the model's population.
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline typename layer::Network<TSeq>::Generator sbm(
    std::vector< size_t > block_sizes,
    std::vector< double > mixing_matrix,
    bool row_major = true
)
{
    return [block_sizes, mixing_matrix, row_major](Model<TSeq> & model) -> AdjList {
        return rgraph_sbm(block_sizes, mixing_matrix, row_major, model);
    };
}

///@}

/**
 * @brief Degree-preserving rewiring, for `layer::Network::set_rewiring()`.
 *
 * @details Each step, a share of the ties, read from the model parameter
 * `param` (a proportion in `[0, 1]`), is rewired: two ties {(i, j), (k, l)}
 * become {(i, l), (k, j)}, skipping swaps that would create self-ties or
 * repeated ties. In a directed network only the sources' ties move, which keeps
 * every agent's in- and out-degree.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline typename layer::Network<TSeq>::Rewiring rewire_degseq(std::string param)
{
    return [ref = ParamRef(std::move(param))](
        layer::Network<TSeq> & net, Model<TSeq> & model
    ) -> void {

        const epiworld_double proportion = ref(model);
        if ((proportion < 0.0) || (proportion > 1.0))
            throw std::range_error(
                "The rewiring proportion '" + ref.name() + "' must be in [0, 1]; "
                "it is " + std::to_string(proportion) + "."
            );

        epiworld::rewire_degseq(net, model, proportion);

    };
}

}

#endif
