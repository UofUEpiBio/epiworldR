#ifndef EPIWORLD_SAMPLER_MIXING_HPP
#define EPIWORLD_SAMPLER_MIXING_HPP

namespace sampler {

/**
 * @brief Samples infectious contacts for models with mixing (contact matrix).
 *
 * @details The sampler is configured with two lists of states: the
 * *infectious* states (whose agents can be contacted and transmit) and the
 * *available* states (agents that take part in the mixing, which is what the
 * contact rate is spread over). Agents are grouped by their first entity.
 *
 * Each step, `update()` rebuilds the pool of infectious agents of each group,
 * and the number of available agents per group, without scanning the
 * population: infectious agents come from the model's per-state index (and are
 * walked in ascending id order), and available counts are the group sizes
 * minus the few agents in non-available states.
 *
 * `sample()` then draws, for each group `g`, a binomial number of contacts
 * (with probability `contact_rate(group of agent, g) / available(g)` over the
 * infectious agents of `g`), picks the infectious agents with replacement
 * (duplicates are kept) and rejects contacts with the agent itself.
 *
 * @tparam TSeq Type of sequence.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class Mixing {
private:

    std::vector< char > is_infectious; ///< Lookup table, one entry per state
    std::vector< char > is_available;  ///< Lookup table, one entry per state
    std::vector< int > infectious_states;
    std::vector< int > available_states;

    std::vector< size_t > pool;      ///< Infectious agents, partitioned by group
    std::vector< size_t > offsets;   ///< Where each group starts in `pool`
    std::vector< size_t > counts;    ///< Infectious agents per group
    std::vector< double > inv_available; ///< 1 / (available agents), per group
    std::vector< size_t > n_available;

    std::vector< uint64_t > marks;   ///< Bitset used to order the pool
    std::vector< size_t > ordered;   ///< Infectious ids (ascending), scratch
    std::vector< size_t > sampled;   ///< Output buffer of `sample()`

public:

    Mixing() = default;

    /**
     * @param infectious States whose agents can be sampled as contacts.
     * @param available States whose agents take part in the mixing.
     */
    Mixing(
        std::vector< int > infectious,
        std::vector< int > available
    ) : infectious_states(std::move(infectious)),
        available_states(std::move(available)) {}

    /**
     * @brief Sizes all buffers for the model and builds the first pools.
     * Call it from the model's `reset()`, after `Model::reset()`.
     */
    void reset(Model<TSeq> & model);

    /**
     * @brief Rebuilds the infectious pools and the available counts from the
     * current state of the model (call it once per step).
     */
    void update(Model<TSeq> & model);

    /**
     * @brief Samples the infectious contacts of `agent`.
     * @return Number of valid entries in `sampled_ids()`.
     */
    size_t sample(
        Agent<TSeq> * agent,
        Model<TSeq> & model,
        const ContactMatrix & contacts
    );

    /// Agent ids drawn by the last `sample()` (first `n` entries are valid).
    const std::vector< size_t > & sampled_ids() const { return sampled; }

    /// Number of infectious agents in a group.
    size_t get_n_infectious(size_t group) const { return counts[group]; }

};

template<typename TSeq>
inline void Mixing<TSeq>::reset(Model<TSeq> & model)
{

    const size_t ns = model.get_n_states();
    is_infectious.assign(ns, 0);
    is_available.assign(ns, 0);

    for (int s : infectious_states)
    {
        if ((s < 0) || (static_cast< size_t >(s) >= ns))
            throw std::range_error("Mixing sampler: infectious state out of range.");
        is_infectious[s] = 1;
    }

    for (int s : available_states)
    {
        if ((s < 0) || (static_cast< size_t >(s) >= ns))
            throw std::range_error("Mixing sampler: available state out of range.");
        is_available[s] = 1;
    }

    const size_t ng = model.get_entities().size();
    const size_t n = model.size();

    pool.assign(n, 0u);
    sampled.assign(n, 0u);
    ordered.clear();
    ordered.reserve(n);
    counts.assign(ng, 0u);
    offsets.assign(ng, 0u);
    n_available.assign(ng, 0u);
    inv_available.assign(ng, 0.0);
    marks.assign((n + 63u) / 64u, 0u);

    // Groups start where the previous one ends (at most its size)
    for (size_t i = 1u; i < ng; ++i)
        offsets[i] = offsets[i - 1u] + model.get_entity(i - 1u).size();

    update(model);

}

template<typename TSeq>
inline void Mixing<TSeq>::update(Model<TSeq> & model)
{

    const size_t ns = is_infectious.size();
    const size_t ng = counts.size();

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

    // Placing them in their group (agents without entity are not in the mixing)
    for (size_t id : ordered)
    {
        auto & a = model.get_agent(id);
        if (a.get_n_entities() == 0u)
            continue;

        const size_t g = a.get_entity(0u, model).get_id();
        pool[offsets[g] + counts[g]++] = id;
    }

    // Available = group size - agents in non-available states
    for (size_t g = 0u; g < ng; ++g)
        n_available[g] = model.get_entity(g).size();

    for (size_t s = 0u; s < ns; ++s)
    {
        if (is_available[s])
            continue;

        for (size_t id : model.get_agents_in_state(s))
        {
            auto & a = model.get_agent(id);
            if (a.get_n_entities() == 0u)
                continue;

            const size_t g = a.get_entity(0u, model).get_id();
            if (n_available[g] > 0u)
                --n_available[g];
        }
    }

    for (size_t g = 0u; g < ng; ++g)
        inv_available[g] = (n_available[g] > 0u) ?
            1.0 / static_cast< double >(n_available[g]) : 0.0;

}

template<typename TSeq>
inline size_t Mixing<TSeq>::sample(
    Agent<TSeq> * agent,
    Model<TSeq> & model,
    const ContactMatrix & contacts
)
{

    const size_t agent_group_id = agent->get_entity(0u, model).get_id();
    const size_t ngroups = counts.size();

    size_t samp_id = 0u;
    for (size_t g = 0u; g < ngroups; ++g)
    {

        const size_t group_size = counts[g];

        if (group_size == 0u)
            continue;

        // How many from this entity?
        int nsamples = model.rbinom(
            group_size,
            inv_available[g] * contacts.get_contact_rate(agent_group_id, g, false)
        );

        for (int s = 0; s < nsamples; ++s)
        {

            // Randomly selecting an agent
            int which = model.runif() * group_size;

            // Correcting overflow error
            if (which >= static_cast< int >(group_size))
                which = static_cast< int >(group_size) - 1;

            const size_t id = pool[offsets[g] + which];

            #ifdef EPI_DEBUG
            if (!is_infectious[model.get_agent(id).get_state()])
                throw std::logic_error(
                    "The agent is not infected, but it should be."
                );
            #endif

            // Can't sample itself
            if (id == static_cast< size_t >(agent->get_id()))
                continue;

            sampled[samp_id++] = id;

        }

    }

    // Reporting the interactions, only if somebody listens
    if (model.has_post_sampling())
        model.register_sampled_contacts(sampled.data(), samp_id, agent->get_id());

    return samp_id;

}

}

#endif
