#ifndef EPIWORLD_LAYER_MIXING_HPP
#define EPIWORLD_LAYER_MIXING_HPP

namespace layer {

template<typename TSeq>
class Mixing;

/**
 * @brief A mixing layer (see `Mixing`).
 * @param name Name of the layer (unique in a model).
 * @ingroup contact_layers
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
std::shared_ptr< Mixing<TSeq> > mixing(std::string name);

/**
 * @brief Contact layer with mixing: contacts drawn fresh every step, by group.
 *
 * @details Created with `layer::mixing()`. Every agent of a group takes part
 * in the mixing (the contacts are spread over all of them), and the agents in
 * the model's infectious states transmit (`Model::set_infectious_states()`).
 * The quarantine model takes its isolated agents out of the mixing altogether
 * (see `Model::set_excluded_states()`).
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
 *   number of contacts per day of an agent in `g` with agents in `h`. With
 *   `set_contact_rate(param)`, the contacts are the model parameter `param`
 *   times the matrix, so a global event or a calibration can change them
 *   (one group and a 1 x 1 matrix of 1: everyone meets "Contact rate" others
 *   a day, at random).
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
 * and the number of agents in the mixing per group, without scanning the
 * population: infectious agents are the carriers in infectious states, from
 * the model's per-state index (walked in ascending id order), and the counts
 * are the group sizes minus the few agents in excluded states.
 *
 * `sample()` then meets each infectious agent of a group `g` with a
 * probability `p` per pair, independently: it draws a binomial number of
 * contacts over the infectious agents of `g` and picks that many *distinct*
 * ones (a partial shuffle of the pool), skipping the agent itself. With
 * `p = 1` it takes them all, with no random draws. Within groups, the only
 * group is the agent's own and `p = q[g]`. With a matrix, `p = min(1,
 * rate * contact_matrix(group of agent, g) / others(g))`, where `rate` is the
 * contact-rate parameter (1 without one) and `others(g)` is the number of
 * agents of `g` in the mixing other than the agent itself: the contacts per
 * day stay as the matrix says when agents are excluded.
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
 * ### History (traced within-group layers)
 *
 * A within-group layer with tracing on (`enable_tracing()`) keeps, for the
 * current run, the history of who was in which group and who was absent (in
 * an excluded state), day by day: `groupmates_of(agent, d1, d2, out)` lists
 * who shared a group with `agent`, both present, between days `d1` and `d2`, with
 * the days they shared. With `q = 1` that is exactly who the agent met; with
 * `q < 1`, who it could have met (contact tracing by membership: "notify the
 * household"). A move or an absence is recorded when it takes effect (the end
 * of the step, so the next day); leaving the layer (`NO_GROUP`) is a move like
 * any other. Agents that never change group or presence add nothing to the
 * history. Editing the grouping between runs drops it.
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

    std::optional< ParamRef > rate_ref; ///< Empty: rate 1

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

    // History of the run (traced within-group layers, see groupmates_of())

    /// From day `from` on, the agent is in `group`, present or absent.
    struct Spell {
        int from;
        bool present;
        size_t group;
    };

    /// An agent left a group on day `day` (its first day elsewhere).
    struct Departure {
        int day;
        size_t agent;
    };

    bool history = false; ///< Kept this run, and still valid
    /// Per agent, its spells in order (sized only with a history)
    std::vector< std::vector< Spell > > spells_of;
    /// Per group, the departures in order (sized only with a history)
    std::vector< std::vector< Departure > > departures;
    std::vector< bool > absent;       ///< Per agent, in an excluded state
    std::vector< size_t > absent_ids; ///< The agents flagged in `absent`
    /// The step's movers and their group before it (scratch of `update()`)
    std::vector< std::pair< size_t, size_t > > moved_from;

    /// Records the agent's group and presence from tomorrow on.
    void record(size_t agent_id, size_t old_group, int today);
    void clear_history();

    /// Calls `f(group, first, last)` for each span of days in `[a, b]` that
    /// the agent spent in a group, present.
    template<typename F>
    void spells(size_t agent_id, int a, int b, F f) const;

    /// The day the agent last left group `g` (it is no longer there).
    int last_departure(size_t agent_id, size_t g) const;

    std::vector< size_t > pool;      ///< Infectious agents, partitioned by group
    std::vector< size_t > offsets;   ///< Where each group starts in `pool`
    std::vector< size_t > counts;    ///< Infectious agents per group
    std::vector< size_t > n_available; ///< Agents in the mixing, per group

    std::vector< uint64_t > marks;   ///< Bitset used to order the pool
    std::vector< size_t > ordered;   ///< Infectious ids (ascending), scratch

    /**
     * @brief How many agents of a group an agent meets in a day: the draw of
     * `Model::rbinom(n, p)` (same distribution, same random numbers), with
     * its parameters set once a day instead of at every draw.
     */
    struct DailyDraw {
        int n = 0;      ///< Infectious agents of the group
        double p = 0.0; ///< Probability of meeting each
        bool poisson = false;
        std::poisson_distribution<>::param_type pois;
        std::binomial_distribution<>::param_type binom;
        void set(int n_, double p_);
        int operator()(Model<TSeq> & model) const;
    };

    /// With a contact matrix: entry `own * n_groups + g` is the draw of an
    /// agent of group `own` (in the mixing) from group `g`, set by `update()`.
    std::vector< DailyDraw > daily;
    void daily_build(const Model<TSeq> & model);

    void rebuild_layout();
    void rebuild_members();
    void move(size_t agent_id, size_t group, bool notify_queue);
    void check_agent(size_t agent_id) const;
    void check_group(size_t group) const;

    /// Meets `n_met` distinct infectious agents of group `g`.
    void draw(size_t g, size_t n_met, size_t agent, Model<TSeq> & model,
        std::vector< size_t > & out);

protected:

    // Created with layer::mixing()
    explicit Mixing(std::string name) : ContactLayer<TSeq>(std::move(name)) {}

    template<typename T>
    friend std::shared_ptr< Mixing<T> > mixing(std::string);

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
     * @brief The contacts are the model parameter `param` times the contact
     * matrix (read at the end of every step, after the global events, so a
     * change to the parameter applies from the next step).
     * @details For a layer with a contact matrix only; `reset()` throws for a
     * layer with contact probabilities.
     */
    void set_contact_rate(std::string param) { rate_ref.emplace(std::move(param)); }

    bool has_contact_rate() const { return rate_ref.has_value(); }

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

    /// A groupmate of an agent over a window of days (see `groupmates_of()`).
    struct Groupmate {
        size_t agent;
        int last_day; ///< The last day they shared a group, both present
        int n_days;   ///< The days they shared a group, both present
    };

    /**
     * @brief Who shared a group with `agent`, both present, between days `d1`
     * and `d2` (traced within-group layers; see the class description).
     * @details One entry per groupmate, however many times either moved or
     * was absent in between. Days after today are not known yet and are left
     * out. The history is that of the current (or last) run, and covers it
     * whole: tracing turned on during a run starts it at the next one.
     * @throws std::logic_error Without a history (see the class description).
     * @param out Replaced by the groupmates, in ascending id order (a vector
     * reused across calls allocates nothing once it has grown).
     */
    void groupmates_of(
        size_t agent_id,
        int d1,
        int d2,
        std::vector< Groupmate > & out
    ) const;

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
     * @details `exclude` is ignored: the pools hold the agents in the model's
     * infectious states only. See `ContactLayer::sample()`.
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
inline std::shared_ptr< Mixing<TSeq> > mixing(std::string name)
{
    return std::shared_ptr< Mixing<TSeq> >(new Mixing<TSeq>(std::move(name)));
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
    clear_history();

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
    clear_history();

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

        if (has_contact_rate())
            throw std::logic_error(
                "Mixing layer '" + this->get_name() + "': a contact rate "
                "scales a contact matrix; this layer has contact "
                "probabilities (within groups)."
            );

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

    clear_history();
    history = within_groups() && this->is_tracing();
    if (history)
    {
        spells_of.resize(n);
        departures.resize(n_groups);
        absent.assign(n, false);
    }

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

    // Moves made during the step; the history keeps each agent's net change
    // (several moves of one agent in a step are one, or none if it is back)
    moved_from.clear();
    for (const auto & [agent_id, group] : moves)
    {
        if (history)
            moved_from.emplace_back(agent_id, groups[agent_id]);

        move(agent_id, group, true);
    }
    moves.clear();

    // The first entry of each agent has its group before the step
    std::stable_sort(
        moved_from.begin(), moved_from.end(),
        [](const auto & x, const auto & y) { return x.first < y.first; }
    );
    for (size_t i = 0u; i < moved_from.size(); ++i)
    {
        const auto & [agent_id, old] = moved_from[i];
        if ((i > 0u) && (moved_from[i - 1u].first == agent_id))
            continue;

        if (groups[agent_id] != old)
            record(agent_id, old, model.today());
    }

    if (layout_stale)
        rebuild_layout();

    const size_t ns = static_cast< size_t >(model.nstates);
    const auto & not_infectious = model.state_not_infectious;

    std::fill(counts.begin(), counts.end(), 0u);

    // Marking the carriers in infectious states; walking the bitset gives
    // ascending ids. States without carriers (e.g., the susceptible) are
    // skipped without a scan.
    for (size_t s = 0u; s < ns; ++s)
    {
        if ((model.state_carriers[s] == 0u) ||
            (!not_infectious.empty() && not_infectious[s]))
            continue;

        for (size_t id : model.state_index_members(s))
            if (model.agent_carrier[id])
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

    // In the mixing = group size - agents in excluded states
    for (size_t g = 0u; g < n_groups; ++g)
        n_available[g] = group_size[g];

    // Back from an excluded state
    if (history)
    {
        size_t k = 0u;
        for (size_t id : absent_ids)
        {
            if (model.state_excluded[model.agent_state[id]])
                absent_ids[k++] = id;
            else
            {
                absent[id] = false;
                record(id, NO_GROUP, model.today());
            }
        }
        absent_ids.resize(k);
    }

    for (size_t s = 0u; s < ns; ++s)
    {
        if (!model.state_excluded[s])
            continue;

        for (size_t id : model.state_index_members(s))
        {
            if (history && !absent[id])
            {
                absent[id] = true;
                absent_ids.push_back(id);
                record(id, NO_GROUP, model.today());
            }

            const size_t g = groups[id];
            if (g == NO_GROUP)
                continue;

            if (n_available[g] > 0u)
                --n_available[g];
        }
    }

    if (!within_groups())
        daily_build(model);

}

template<typename TSeq>
inline void Mixing<TSeq>::record(size_t agent_id, size_t old_group, int today)
{

    const size_t g = groups[agent_id];
    spells_of[agent_id].push_back({today + 1, !absent[agent_id], g});

    // Absences do not change the members
    if ((old_group != NO_GROUP) && (old_group != g))
        departures[old_group].push_back({today + 1, agent_id});

}

template<typename TSeq>
inline void Mixing<TSeq>::clear_history()
{
    history = false;
    // Freed, not just emptied: an untraced layer holds nothing
    std::vector< std::vector< Spell > >().swap(spells_of);
    std::vector< std::vector< Departure > >().swap(departures);
    std::vector< bool >().swap(absent);
    absent_ids.clear();
}

template<typename TSeq>
template<typename F>
inline void Mixing<TSeq>::spells(size_t agent_id, int a, int b, F f) const
{

    // Walking back from the last spell, each one holding until the next
    int end = b + 1;
    const auto & sp = spells_of[agent_id];
    for (auto s = sp.rbegin(); (s != sp.rend()) && (end > a); ++s)
    {
        const int first = std::max(s->from, a);
        if (s->present && (s->group != NO_GROUP) && (first < end))
            f(s->group, first, end - 1);

        end = std::min(end, s->from);
    }

    // Before its first change, the agent is where the run started
    if ((end > a) && (groups_baseline[agent_id] != NO_GROUP))
        f(groups_baseline[agent_id], a, end - 1);

}

template<typename TSeq>
inline int Mixing<TSeq>::last_departure(size_t agent_id, size_t g) const
{

    // The spell after its last one in g (none: it left its starting group)
    const auto & sp = spells_of[agent_id];
    for (size_t k = sp.size() - 1u; k > 0u; --k)
        if (sp[k - 1u].group == g)
            return sp[k].from;

    return sp[0u].from;

}

template<typename TSeq>
inline void Mixing<TSeq>::groupmates_of(
    size_t agent_id,
    int d1,
    int d2,
    std::vector< Groupmate > & out
) const
{

    check_agent(agent_id);
    if (!history)
        throw std::logic_error(
            "Mixing layer '" + this->get_name() + "': " + (
                (this->is_tracing() && within_groups()) ?
                "no history. It covers whole runs: tracing was turned on, "
                "or the grouping edited, after the last run started; the "
                "next run keeps one." :
                "groupmates_of() needs a traced layer (enable_tracing() "
                "before the run) with contact probabilities (within groups)."
            )
        );

    d1 = std::max(d1, 0);
    d2 = std::min(d2, this->model->today());

    // One entry per shared span, merged per groupmate below
    out.clear();
    spells(agent_id, d1, d2, [&](size_t g, int a, int b) {

        auto shared = [&](size_t other) {
            spells(other, a, b, [&](size_t h, int first, int last) {
                if (h == g)
                    out.push_back({other, last, last - first + 1});
            });
        };

        // The members now
        for (size_t other : members[g])
            if (other != agent_id)
                shared(other);

        // Those who left after day a, each once: at its last departure
        const auto & left = departures[g];
        for (auto d = left.rbegin(); (d != left.rend()) && (d->day > a); ++d)
            if (
                (d->agent != agent_id) && (groups[d->agent] != g) &&
                (last_departure(d->agent, g) == d->day)
            )
                shared(d->agent);

    });

    std::sort(out.begin(), out.end(), [](const auto & x, const auto & y) {
        return x.agent < y.agent;
    });

    size_t k = 0u;
    for (size_t i = 0u; i < out.size(); ++i)
    {
        if ((k > 0u) && (out[k - 1u].agent == out[i].agent))
        {
            out[k - 1u].last_day = std::max(out[k - 1u].last_day, out[i].last_day);
            out[k - 1u].n_days  += out[i].n_days;
        }
        else
            out[k++] = out[i];
    }
    out.resize(k);

}

template<typename TSeq>
inline void Mixing<TSeq>::DailyDraw::set(int n_, double p_)
{

    n = n_;
    p = p_;
    poisson = false;

    if ((n == 0) || (p <= 0.0) || (p >= 1.0))
        return;

    // The same choice as Model::rbinom(n, p)
    #ifdef EPI_FAST_BINOM
    if ((p <= 0.01) && (static_cast< double >(n) * p * p <= 0.1))
    {
        poisson = true;
        pois = std::poisson_distribution<>::param_type(static_cast< double >(n) * p);
        return;
    }
    #endif

    binom = std::binomial_distribution<>::param_type(n, p);

}

template<typename TSeq>
inline int Mixing<TSeq>::DailyDraw::operator()(Model<TSeq> & model) const
{

    if ((n == 0) || (p <= 0.0))
        return 0;

    // Everyone, with no random draws
    if (p >= 1.0)
        return n;

    if (poisson)
        return std::min(model.rpoissd(*model.engine, pois), n);

    return model.rbinomd(*model.engine, binom);

}

template<typename TSeq>
inline void Mixing<TSeq>::daily_build(const Model<TSeq> & model)
{

    // The contact rate scales the whole matrix
    const double scale = rate_ref ? (*rate_ref)(model) : 1.0;

    daily.resize(n_groups * n_groups);
    for (size_t own = 0u; own < n_groups; ++own)
    {

        // Column of the group: entry (own, g) is rate[g * n_groups]
        const double * rate = contact_matrix.data() + own;
        for (size_t g = 0u; g < n_groups; ++g)
        {

            // In its own group, an agent's contacts are the others
            const size_t others = n_available[g] - ((g == own) ? 1u : 0u);
            const double p = (others == 0u || n_available[g] == 0u) ? 0.0 :
                std::min(
                    1.0,
                    scale * rate[g * n_groups] / static_cast< double >(others)
                );

            daily[own * n_groups + g].set(static_cast< int >(counts[g]), p);

        }

    }

}

template<typename TSeq>
inline void Mixing<TSeq>::draw(
    size_t g,
    size_t n_met,
    size_t agent,
    Model<TSeq> & model,
    std::vector< size_t > & out
)
{

    const size_t n_pool = counts[g];

    // A binomial number of distinct agents is the same as meeting each one
    // independently with probability p; with p = 1 (a clique) it is everyone,
    // with no draws. A partial shuffle of the pool puts the ones met first.
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
        if (model.get_agent(id).get_virus() == nullptr)
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
        const double q = contact_prob[own];
        const size_t n_pool = counts[own];
        if ((n_pool > 0u) && (q > 0.0))
            draw(
                own,
                (q < 1.0) ?
                    static_cast< size_t >(model.rbinom(static_cast< int >(n_pool), q)) :
                    n_pool,
                agent, model, out
            );
    }
    else if (!model.state_excluded[model.agent_state[agent]])
    {

        // The day's draws of the agent's group (see update())
        const DailyDraw * day = daily.data() + own * n_groups;
        for (size_t g = 0u; g < n_groups; ++g)
        {
            const int n_met = day[g](model);
            if (n_met > 0)
                draw(g, static_cast< size_t >(n_met), agent, model, out);
        }

    }
    else
    {

        // An excluded agent is not one of the others of its group (it rarely
        // samples: its state usually has no sampler)
        const double scale = rate_ref ? (*rate_ref)(model) : 1.0;
        const double * rate = contact_matrix.data() + own;
        for (size_t g = 0u; g < n_groups; ++g)
        {
            const size_t n_pool = counts[g];
            const double p = (n_available[g] == 0u) ? 0.0 : std::min(
                1.0,
                scale * rate[g * n_groups] / static_cast< double >(n_available[g])
            );
            if ((n_pool == 0u) || (p <= 0.0))
                continue;

            draw(
                g,
                (p < 1.0) ?
                    static_cast< size_t >(model.rbinom(static_cast< int >(n_pool), p)) :
                    n_pool,
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
