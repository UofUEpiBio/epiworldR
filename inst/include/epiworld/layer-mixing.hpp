#ifndef EPIWORLD_LAYER_MIXING_HPP
#define EPIWORLD_LAYER_MIXING_HPP

namespace layer {

template<typename TSeq>
class Mixing;

/**
 * @brief A mixing layer (see `Mixing`).
 * @param name Name of the layer (unique in a model).
 * @param infectious States whose agents can be sampled as contacts.
 * @param available States whose agents take part in the mixing.
 * @ingroup contact_layers
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
std::shared_ptr< Mixing<TSeq> > mixing(
    std::string name,
    std::vector< int > infectious,
    std::vector< int > available
);

/**
 * @brief Contact layer with mixing: contacts drawn fresh every step, by group.
 *
 * @details Created with `layer::mixing()`. The layer is configured with two
 * lists of states: the
 * *infectious* states (whose agents can be contacted and transmit) and the
 * *available* states (agents that take part in the mixing, which is what the
 * contact rate is spread over).
 *
 * ### Groups
 *
 * The layer owns the grouping: `set_groups()` assigns every agent a group (or
 * `NO_GROUP`, which leaves it out of the layer). If no grouping is set, the
 * layer builds one from the model's entities once they are distributed (at the
 * first `update()` of a run), each agent going to the group of its first entity
 * (and agents without entities to none), which is how the mixing models are
 * built. Either way, from then on the layer works on its own table, not on the
 * entities.
 *
 * ### Contacts: a matrix or a probability per pair
 *
 * The contacts are given in one of two ways:
 * - `set_contact_matrix()`: contacts across groups, for a few groups that meet
 *   each other (age groups in public spaces). Entry `(g, h)` is the expected
 *   number of contacts per day of an agent in `g` with agents in `h`.
 * - `set_contact_prob()`: contacts within groups only, for many small groups
 *   (households, classrooms, transit trips). Each member meets each other
 *   member of its group with probability `q` per day, whatever the size of the
 *   group; `q = 1` makes every group a clique (everyone meets everyone, every
 *   day), with no ties stored. One number sets the same `q` for every group.
 *   It needs the grouping from `set_groups()`.
 *
 * Setting one between runs drops the other. A within-group layer is seen by
 * the model's queue (an infectious agent reaches only its own group), and
 * `members_of()` lists a group.
 *
 * ### Sampling
 *
 * Each step, `update()` rebuilds the pool of infectious agents of each group,
 * and the number of available agents per group, without scanning the
 * population: infectious agents come from the model's per-state index (and are
 * walked in ascending id order), and available counts are the group sizes
 * minus the few agents in non-available states.
 *
 * `sample()` then meets each infectious agent of a group `g` with a
 * probability `p` per pair, independently: it draws a binomial number of
 * contacts over the infectious agents of `g` and picks that many *distinct*
 * ones (a partial shuffle of the pool), skipping the agent itself. With
 * `p = 1` it takes them all, with no random draws. Within groups, the only
 * group is the agent's own and `p = q[g]`. With a matrix, `p = min(1,
 * contact_matrix(group of agent, g) / others(g))`, where `others(g)` is the
 * number of available agents of `g` other than the agent itself: the
 * contacts per day stay as the matrix says when agents become unavailable.
 *
 * ### Edits
 *
 * `set_group()`, `set_contact_matrix()` and `set_contact_prob()` called
 * during a run apply to that run only: they set the dirty flag, and the next
 * `reset()` restores the baseline. A move made during a step takes effect at
 * the end of the step, with the pools (the agents of the step sample from the
 * grouping it started with). Called between runs, they set the baseline (so a
 * calibration can reuse one model). `set_groups()` sets the baseline grouping
 * and can only be called between runs.
 *
 * @tparam TSeq Type of sequence.
 * @ingroup contact_layers
 */
template<typename TSeq>
class Mixing : public ContactLayer<TSeq>
{
public:

    /// Group of the agents that are not in the layer.
    static constexpr size_t NO_GROUP = std::numeric_limits< size_t >::max();

private:

    std::vector< char > is_infectious; ///< Lookup table, one entry per state
    std::vector< char > is_available;  ///< Lookup table, one entry per state
    std::vector< int > infectious_states;
    std::vector< int > available_states;

    // Grouping
    std::vector< size_t > groups;          ///< Group of each agent
    std::vector< size_t > groups_baseline; ///< Set by set_groups()
    size_t n_groups = 0u;
    bool explicit_groups = false; ///< Set by set_groups(), else from entities
    bool groups_pending = false;  ///< Groups to read from the entities
    bool layout_stale = true;     ///< Group sizes and pool offsets to rebuild
    std::vector< size_t > group_size; ///< Agents in each group

    std::vector< double > contact_matrix;          ///< Column-major
    std::vector< double > contact_matrix_baseline; ///< Set between runs
    std::vector< double > contact_prob;          ///< Within-group `q`, per group
    std::vector< double > contact_prob_baseline; ///< Set between runs (one entry: every group)

    std::vector< std::vector< size_t > > members; ///< Agents of each group

    /// Moves made during a step (agent, group), applied by `update()`
    std::vector< std::pair< size_t, size_t > > moves;

    std::vector< size_t > pool;      ///< Infectious agents, partitioned by group
    std::vector< size_t > offsets;   ///< Where each group starts in `pool`
    std::vector< size_t > counts;    ///< Infectious agents per group
    std::vector< size_t > n_available; ///< Available agents, per group

    std::vector< uint64_t > marks;   ///< Bitset used to order the pool
    std::vector< size_t > ordered;   ///< Infectious ids (ascending), scratch

    void rebuild_layout();
    void rebuild_members();
    void move(size_t agent_id, size_t group, bool notify_queue);
    void check_agent(size_t agent_id) const;
    void check_group(size_t group) const;

    /// Meets each infectious agent of group `g` with probability `p`.
    void draw(size_t g, double p, size_t agent, Model<TSeq> & model,
        std::vector< size_t > & out);

protected:

    // Created with layer::mixing()
    Mixing(
        std::string name,
        std::vector< int > infectious,
        std::vector< int > available
    ) : ContactLayer<TSeq>(std::move(name)),
        infectious_states(std::move(infectious)),
        available_states(std::move(available)) {}

    template<typename T>
    friend std::shared_ptr< Mixing<T> > mixing(
        std::string, std::vector< int >, std::vector< int >
    );

    // Reads the groups unchecked (as the ties of a network layer)
    friend class Queue<TSeq>;

public:

    /**
     * @brief Sets the grouping of the layer (between runs only).
     * @param group_of_agent One entry per agent: a group in `[0, n_groups)`, or
     * `NO_GROUP` to leave the agent out of the layer.
     * @param n_groups_ Number of groups (the contact matrix is
     * `n_groups x n_groups`).
     */
    void set_groups(std::vector< size_t > group_of_agent, size_t n_groups_);

    /// Group of an agent (`NO_GROUP` if it is not in the layer).
    size_t get_group(size_t agent_id) const
    {
        check_agent(agent_id);
        return groups[agent_id];
    }

    /// Whether the agent takes part in the layer (valid after `reset()`).
    bool in_layer(size_t agent_id) const { return groups[agent_id] != NO_GROUP; }

    /**
     * @brief Moves an agent to another group (`NO_GROUP` takes it out of the
     * layer).
     * @details Like every change made during a step, it is in force for the
     * contacts of the next step: the move waits for the end of the step, when
     * the pools are rebuilt (until then, `get_group()` gives the old group).
     */
    void set_group(size_t agent_id, size_t group);

    size_t get_n_groups() const { return n_groups; }

    /**
     * @brief Sets the contact matrix.
     * @param cmat Column-major `n_groups x n_groups` matrix: entry `(i, j)`
     * is the expected number of contacts per day that an agent in group `i` has
     * with agents in group `j`.
     */
    void set_contact_matrix(std::vector< double > cmat);

    const std::vector< double > & get_contact_matrix() const
    {
        return contact_matrix;
    }

    /**
     * @brief Sets the probability that two members of a group meet on a given
     * day (contacts within groups only).
     * @param q One probability per group, or a single one for every group;
     * `1` makes the groups cliques (everyone meets everyone, every day).
     */
    void set_contact_prob(std::vector< double > q);
    void set_contact_prob(double q)
    {
        set_contact_prob(std::vector< double >{q});
    }

    /// The within-group contact probabilities (one per group after `reset()`;
    /// empty with a contact matrix).
    const std::vector< double > & get_contact_prob() const
    {
        return contact_prob;
    }

    /// Whether the layer mixes within groups only (`set_contact_prob()`).
    bool within_groups() const { return !contact_prob_baseline.empty(); }

    /// The agents of group `g`, in no particular order (valid after `reset()`).
    const std::vector< size_t > & members_of(size_t g) const
    {
        check_group(g);
        return members[g];
    }

    layer::Mixing<TSeq> * as_mixing() override { return this; }

    /// Within groups only: a carrier reaches only its own group.
    bool supports_queue() const override { return within_groups(); }

    /**
     * @brief Restores the baseline and sizes all buffers for the model (the
     * first `update()`, on day 0, builds the pools).
     */
    void reset(Model<TSeq> & model) override;

    /**
     * @brief Rebuilds the infectious pools and the available counts from the
     * current state of the model (once per step).
     */
    void update(Model<TSeq> & model) override;

    std::shared_ptr< ContactLayer<TSeq> > clone() const override
    {
        return std::make_shared< Mixing<TSeq> >(*this);
    }

    /**
     * @brief Samples the infectious contacts of `agent` (none if it is not in
     * the layer), appending them to `out`.
     * @details `exclude` is ignored: the infectious states of the layer say
     * who can be contacted. See `ContactLayer::sample()`.
     */
    size_t sample(
        size_t agent,
        Model<TSeq> & model,
        const std::vector< bool > * exclude,
        std::vector< size_t > & out
    ) override;

    /// Number of infectious agents in a group.
    size_t get_n_infectious(size_t group) const { return counts[group]; }

};

template<typename TSeq>
inline std::shared_ptr< Mixing<TSeq> > mixing(
    std::string name,
    std::vector< int > infectious,
    std::vector< int > available
)
{
    return std::shared_ptr< Mixing<TSeq> >(new Mixing<TSeq>(
        std::move(name), std::move(infectious), std::move(available)
    ));
}

template<typename TSeq>
inline void Mixing<TSeq>::set_groups(
    std::vector< size_t > group_of_agent,
    size_t n_groups_
)
{

    for (size_t g : group_of_agent)
        if ((g != NO_GROUP) && (g >= n_groups_))
            throw std::range_error(
                "Mixing layer '" + this->get_name() + "': group " +
                std::to_string(g) + " is out of range (" +
                std::to_string(n_groups_) + " groups)."
            );

    if (!this->edits_baseline())
        throw std::logic_error(
            "Mixing layer '" + this->get_name() + "': set_groups() sets the "
            "baseline and cannot be called during a run; use set_group()."
        );

    groups_baseline = std::move(group_of_agent);
    groups          = groups_baseline;
    n_groups        = n_groups_;
    explicit_groups = true;
    layout_stale    = true;
    rebuild_members();

}

template<typename TSeq>
inline void Mixing<TSeq>::check_agent(size_t agent_id) const
{
    if (agent_id >= groups.size())
        throw std::range_error(
            "Mixing layer '" + this->get_name() + "': agent " +
            std::to_string(agent_id) + " is out of range (the layer has " +
            std::to_string(groups.size()) + " agents; was it reset?)."
        );
}

template<typename TSeq>
inline void Mixing<TSeq>::check_group(size_t group) const
{
    if (group >= n_groups)
        throw std::range_error(
            "Mixing layer '" + this->get_name() + "': group " +
            std::to_string(group) + " is out of range (" +
            std::to_string(n_groups) + " groups)."
        );
}

template<typename TSeq>
inline void Mixing<TSeq>::set_group(size_t agent_id, size_t group)
{

    check_agent(agent_id);
    if (group != NO_GROUP)
        check_group(group);

    // During a run, the move waits for the end of the step (see update())
    if (!this->edits_baseline())
    {
        moves.emplace_back(agent_id, group);
        return;
    }

    if (!explicit_groups)
        throw std::logic_error(
            "Mixing layer '" + this->get_name() + "': set_groups() must be "
            "called before set_group() outside a run."
        );

    groups_baseline[agent_id] = group;
    move(agent_id, group, false);

}

template<typename TSeq>
inline void Mixing<TSeq>::move(
    size_t agent_id,
    size_t group,
    bool notify_queue
)
{

    const size_t old = groups[agent_id];
    if (old == group)
        return;

    // During a run, the queue follows the move (it sees within-group layers)
    if (notify_queue && (this->model != nullptr) && within_groups())
        this->model->group_member_moved(
            agent_id,
            old   == NO_GROUP ? nullptr : &members[old],
            group == NO_GROUP ? nullptr : &members[group]
        );

    if (old != NO_GROUP)
    {
        auto & m = members[old];
        auto it = std::find(m.begin(), m.end(), agent_id);
        *it = m.back();
        m.pop_back();
    }

    if (group != NO_GROUP)
        members[group].push_back(agent_id);

    groups[agent_id] = group;
    layout_stale = true;

}

template<typename TSeq>
inline void Mixing<TSeq>::set_contact_matrix(std::vector< double > cmat)
{

    for (double c : cmat)
        if (c < 0.0)
            throw std::range_error(
                "Mixing layer '" + this->get_name() + "': the contact matrix "
                "must be non-negative (" + std::to_string(c) + " < 0)."
            );

    // Between runs, the size is checked by reset() (the grouping may not be
    // set yet); during a run, the grouping and the form are fixed
    if (this->edits_baseline())
    {
        contact_matrix_baseline = cmat;
        contact_prob_baseline.clear();
    }
    else if (within_groups())
        throw std::logic_error(
            "Mixing layer '" + this->get_name() + "': the layer mixes within "
            "groups (set_contact_prob()); a contact matrix can only replace "
            "the contact probabilities between runs."
        );
    else if (cmat.size() != n_groups * n_groups)
        throw std::length_error(
            "Mixing layer '" + this->get_name() + "': the contact matrix has " +
            std::to_string(cmat.size()) + " entries, but the layer has " +
            std::to_string(n_groups) + " groups."
        );

    contact_matrix = std::move(cmat);

}

template<typename TSeq>
inline void Mixing<TSeq>::set_contact_prob(std::vector< double > q)
{

    if (q.empty())
        throw std::length_error(
            "Mixing layer '" + this->get_name() + "': set_contact_prob() needs "
            "at least one probability."
        );

    for (double x : q)
        if (!((x >= 0.0) && (x <= 1.0)))
            throw std::range_error(
                "Mixing layer '" + this->get_name() + "': contact probabilities "
                "must be in [0, 1] (" + std::to_string(x) + ")."
            );

    // Between runs, the size is checked by reset(); during a run, the
    // grouping and the form are fixed
    if (this->edits_baseline())
    {
        contact_prob_baseline = q;
        contact_matrix_baseline.clear();
    }
    else
    {

        if (!within_groups())
            throw std::logic_error(
                "Mixing layer '" + this->get_name() + "': the layer has a "
                "contact matrix; contact probabilities can only replace it "
                "between runs."
            );

        // A copy: assign() may free the element before reading it
        if (q.size() == 1u)
            q.assign(n_groups, double(q[0u]));
        else if (q.size() != n_groups)
            throw std::length_error(
                "Mixing layer '" + this->get_name() + "': " +
                std::to_string(q.size()) + " contact probabilities, but the "
                "layer has " + std::to_string(n_groups) + " groups."
            );

    }

    contact_prob = std::move(q);

}

template<typename TSeq>
inline void Mixing<TSeq>::rebuild_members()
{

    members.assign(n_groups, {});
    for (size_t i = 0u; i < groups.size(); ++i)
        if (groups[i] != NO_GROUP)
            members[groups[i]].push_back(i);

}

template<typename TSeq>
inline void Mixing<TSeq>::rebuild_layout()
{

    group_size.assign(n_groups, 0u);
    for (size_t g : groups)
        if (g != NO_GROUP)
            ++group_size[g];

    counts.assign(n_groups, 0u);
    offsets.assign(n_groups, 0u);
    n_available.assign(n_groups, 0u);

    // Groups start where the previous one ends (at most its size)
    for (size_t i = 1u; i < n_groups; ++i)
        offsets[i] = offsets[i - 1u] + group_size[i - 1u];

    layout_stale = false;

}

template<typename TSeq>
inline void Mixing<TSeq>::reset(Model<TSeq> & model)
{

    const size_t ns = model.get_n_states();
    is_infectious.assign(ns, 0);
    is_available.assign(ns, 0);

    for (int s : infectious_states)
    {
        if ((s < 0) || (static_cast< size_t >(s) >= ns))
            throw std::range_error("Mixing layer: infectious state out of range.");
        is_infectious[s] = 1;
    }

    for (int s : available_states)
    {
        if ((s < 0) || (static_cast< size_t >(s) >= ns))
            throw std::range_error("Mixing layer: available state out of range.");
        is_available[s] = 1;
    }

    const size_t n = model.size();

    if (explicit_groups)
    {

        // Edits of the last run are dropped
        groups = groups_baseline;

        if (groups.size() != n)
            throw std::length_error(
                "Mixing layer '" + this->get_name() + "': the grouping has " +
                std::to_string(groups.size()) + " agents, but the model has " +
                std::to_string(n) + "."
            );

    }
    else
    {

        // One group per entity; the agents are placed at the first update,
        // once the entities are distributed
        n_groups = model.get_entities().size();
        groups.assign(n, NO_GROUP);
        groups_pending = true;

    }

    // Edits of the last run are dropped
    if (within_groups())
    {

        // The queue registers carriers before the first update, so the
        // members must be known now, not taken from the entities later
        if (!explicit_groups)
            throw std::logic_error(
                "Mixing layer '" + this->get_name() + "': a layer with contact "
                "probabilities (within groups) needs its groups from "
                "set_groups()."
            );

        contact_prob = contact_prob_baseline;
        // A copy: assign() may free the element before reading it
        if (contact_prob.size() == 1u)
            contact_prob.assign(n_groups, double(contact_prob[0u]));
        else if (contact_prob.size() != n_groups)
            throw std::length_error(
                "Mixing layer '" + this->get_name() + "': " +
                std::to_string(contact_prob.size()) + " contact "
                "probabilities, but the layer has " +
                std::to_string(n_groups) + " groups."
            );

        contact_matrix.clear();

    }
    else
    {

        contact_matrix = contact_matrix_baseline;
        if (contact_matrix.size() != n_groups * n_groups)
            throw std::length_error(
                "Mixing layer '" + this->get_name() + "': the contact matrix "
                "has " + std::to_string(contact_matrix.size()) + " entries, but "
                "the layer has " + std::to_string(n_groups) + " groups."
            );

        contact_prob.clear();

    }

    rebuild_members();
    moves.clear();

    pool.assign(n, 0u);
    ordered.clear();
    ordered.reserve(n);
    marks.assign((n + 63u) / 64u, 0u);

    rebuild_layout();
    this->mark_clean();

}

template<typename TSeq>
inline void Mixing<TSeq>::update(Model<TSeq> & model)
{

    // The agent goes to the group of its first entity
    if (groups_pending)
    {
        for (size_t i = 0u; i < groups.size(); ++i)
        {
            auto & a = model.get_agent(i);
            if (a.get_n_entities() != 0u)
                groups[i] = a.get_entity(0u, model).get_id();
        }

        groups_pending = false;
        layout_stale = true;
        rebuild_members();
    }

    // Moves made during the step
    for (const auto & [agent_id, group] : moves)
        move(agent_id, group, true);
    moves.clear();

    if (layout_stale)
        rebuild_layout();

    const size_t ns = is_infectious.size();

    std::fill(counts.begin(), counts.end(), 0u);

    // Marking the infectious agents; walking the bitset gives ascending ids
    for (size_t s = 0u; s < ns; ++s)
    {
        if (!is_infectious[s])
            continue;

        for (size_t id : model.get_agents_in_state(s))
            marks[id >> 6] |= (uint64_t(1) << (id & 63u));
    }

    ordered.clear();
    for (size_t w = 0u; w < marks.size(); ++w)
    {
        uint64_t word = marks[w];
        marks[w] = 0u;
        while (word != 0u)
        {
            ordered.push_back((w << 6) + epi_ctz64(word));
            word &= (word - 1u);
        }
    }

    // Placing them in their group (agents out of the layer are not in the mixing)
    for (size_t id : ordered)
    {
        const size_t g = groups[id];
        if (g == NO_GROUP)
            continue;

        pool[offsets[g] + counts[g]++] = id;
    }

    // Available = group size - agents in non-available states
    for (size_t g = 0u; g < n_groups; ++g)
        n_available[g] = group_size[g];

    for (size_t s = 0u; s < ns; ++s)
    {
        if (is_available[s])
            continue;

        for (size_t id : model.get_agents_in_state(s))
        {
            const size_t g = groups[id];
            if (g == NO_GROUP)
                continue;

            if (n_available[g] > 0u)
                --n_available[g];
        }
    }

}

template<typename TSeq>
inline void Mixing<TSeq>::draw(
    size_t g,
    double p,
    size_t agent,
    Model<TSeq> & model,
    std::vector< size_t > & out
)
{

    const size_t n_pool = counts[g];

    if ((n_pool == 0u) || (p <= 0.0))
        return;

    // A binomial number of distinct agents is the same as meeting each one
    // independently with probability p; with p = 1 (a clique) it is everyone,
    // with no draws. A partial shuffle of the pool puts the ones met first.
    const size_t n_met = (p < 1.0) ?
        static_cast< size_t >(model.rbinom(static_cast< int >(n_pool), p)) :
        n_pool;

    size_t * first = pool.data() + offsets[g];
    for (size_t s = 0u; s < n_met; ++s)
    {

        if (n_met < n_pool)
        {
            size_t which = s + static_cast< size_t >(
                model.runif() * static_cast< double >(n_pool - s)
            );

            // Correcting overflow error
            if (which >= n_pool)
                which = n_pool - 1u;

            std::swap(first[s], first[which]);
        }

        const size_t id = first[s];

        #ifdef EPI_DEBUG
        if (!is_infectious[model.get_agent(id).get_state()])
            throw std::logic_error(
                "The agent is not infected, but it should be."
            );
        #endif

        // Can't meet itself
        if (id != agent)
            out.push_back(id);

    }

}

template<typename TSeq>
inline size_t Mixing<TSeq>::sample(
    size_t agent,
    Model<TSeq> & model,
    const std::vector< bool > *,
    std::vector< size_t > & out
)
{

    if (!in_layer(agent))
        return 0u;

    const size_t first = out.size();

    const size_t own = groups[agent];

    if (within_groups())
    {
        // Only the agent's own group
        draw(own, contact_prob[own], agent, model, out);
    }
    else
    {

        // In its own group, an available agent's contacts are the others
        const bool available =
            is_available[model.get_agent(agent).get_state()];

        // Column of the agent's group: entry (agent's group, g) is rate[g]
        const double * rate = contact_matrix.data() + own;
        for (size_t g = 0u; g < n_groups; ++g)
        {
            const size_t others =
                n_available[g] - (((g == own) && available) ? 1u : 0u);

            if (others > 0u)
                draw(
                    g,
                    std::min(1.0, rate[g * n_groups] / static_cast< double >(others)),
                    agent, model, out
                );
        }

    }

    // Reporting the interactions, only if somebody listens
    const size_t n_sampled = out.size() - first;
    this->trace(out.data() + first, n_sampled, agent);
    if (model.has_post_sampling())
        model.register_sampled_contacts(
            out.data() + first, n_sampled, static_cast< int >(agent)
        );

    return n_sampled;

}

}

#endif
