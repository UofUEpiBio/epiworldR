#ifndef EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP
#define EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP

/**
 * @file model-meat-transmission.hpp
 * @brief Network transmission by pushing infection odds.
 *
 * @details A susceptible agent `i` whose update function is
 * `default_update_susceptible` (or `sampler::make_update_susceptible()`)
 * *pulls*: it scans its neighbors, collects the per-contact probabilities
 * `p_ij` of those carrying a virus, and `roulette()` draws "no infection" or a
 * single infector. That draw depends only on the odds `r_ij = p_ij / (1 - p_ij)`:
 *
 *     P(no infection) = 1 / (1 + R_i),   P(infected by j) = r_ij / (1 + R_i),
 *
 * with `R_i` the sum of the odds. So the same outcome can be *pushed*: every
 * agent carrying a virus adds its odds to its susceptible neighbors, each of
 * which keeps a weighted reservoir sample of its infectors; then each touched
 * agent makes one draw. The cost follows the carriers' ties instead of the
 * susceptibles', which is much cheaper while an outbreak is small.
 *
 * Contacts with `p_ij >= 1` are certain: the infection happens and the
 * infector is drawn uniformly among them (the limit of the odds as `p -> 1`).
 *
 * The docs page "Push and Pull Transmission" (docs/impl/transmission-sampling.md)
 * has the derivation.
 */

template<typename TSeq>
inline bool Model<TSeq>::transmission_prepare()
{

    const size_t ns = static_cast< size_t >(nstates);

    push_pushable.assign(ns, 0);
    push_default.assign(ns, 0);
    push_excluded.assign(ns * ns, 0);
    push_source_ok.assign(ns, 0);

    bool any = false;

    // Recognizing the samplers needs std::function::target(), i.e., RTTI.
    // Without it, every state keeps pulling.
    #if defined(__cpp_rtti) || defined(__GXX_RTTI) || defined(_CPPRTTI)
    for (size_t t = 0u; t < ns; ++t)
    {

        const auto & fun = state_fun[t];
        if (!fun)
            continue;

        using FunPtr = void(*)(Agent<TSeq>*, Model<TSeq>*);
        const FunPtr * fp = fun.template target< FunPtr >();

        if ((fp != nullptr) && (*fp == &default_update_susceptible<TSeq>))
        {
            push_pushable[t] = 1;
            push_default[t] = 1;
            any = true;
        }
        else if (
            const auto * us = fun.template target< sampler::UpdateSusceptible<TSeq> >()
        )
        {

            push_pushable[t] = 1;
            any = true;

            for (auto s : us->exclude)
            {
                if (s >= ns)
                    throw std::logic_error(
                        std::string("You are trying to exclude a state that is out of range: ") +
                        std::to_string(s) + std::string(". There are only ") +
                        std::to_string(ns) + std::string(" states in the model.")
                        );

                push_excluded[t * ns + s] = 1;
            }

        }

    }
    #endif

    // A state can be a source if some pushable state takes infections from it
    for (size_t t = 0u; t < ns; ++t)
        if (push_pushable[t])
            for (size_t s = 0u; s < ns; ++s)
                if (!push_excluded[t * ns + s])
                    push_source_ok[s] = 1;

    return any;

}

template<typename TSeq>
inline bool Model<TSeq>::transmission_choose_push() const
{

    if (transmission_mode == TransmissionMode::pull)
        return false;

    if (transmission_mode == TransmissionMode::push)
        return true;

    // Pushing walks every tie of every carrier that can transmit; pulling walks
    // every tie of every susceptible agent. Both sums are kept per state, so
    // this is O(number of states). It deliberately ignores the queue: the
    // decision -- and so the random stream -- is the same with queuing on or
    // off. The queue does make pulling cheaper than this sum suggests (it
    // skips susceptibles with no infectious neighbor), which is what kappa < 1
    // accounts for.
    double cost_push = 0.0;
    double cost_pull = 0.0;
    for (size_t s = 0u; s < static_cast< size_t >(nstates); ++s)
    {
        if (push_source_ok[s])
            cost_push += static_cast< double >(state_carrier_degree[s]);
        if (push_pushable[s])
            cost_pull += static_cast< double >(state_degree[s]);
    }

    return cost_push <= transmission_kappa * cost_pull;

}

template<typename TSeq>
inline void Model<TSeq>::transmission_push()
{

    const size_t ns = static_cast< size_t >(nstates);

    if (push_slot.size() != population.size())
        push_slot.assign(population.size(), -1);

    push_targets.clear();

    // Phase 1: every carrier that can transmit adds its odds to its eligible
    // neighbors. Nothing changes state until events_run(), so this sees the
    // model as it was at the start of the step, as pulling does.
    //
    // The carriers are visited in ascending id order (marked in a bitset, then
    // walked), not in the index's order: networks are usually built with
    // neighbors close in id, so this sweeps memory the way a pull does, and on
    // large populations it saves most of the cache misses.
    const size_t nwords = (population.size() + 63u) / 64u;
    if (push_sources.size() != nwords)
        push_sources.assign(nwords, 0u);

    for (size_t s = 0u; s < ns; ++s)
    {

        if (!push_source_ok[s] || (state_carriers[s] == 0u))
            continue;

        for (size_t id : state_index_members(s))
            push_sources[id >> 6] |= (uint64_t(1) << (id & 63u));

    }

    for (size_t w = 0u; w < nwords; ++w)
    {

        uint64_t word = push_sources[w];
        push_sources[w] = 0u;

        while (word != 0u)
        {

            const size_t j_id = (w << 6) + epi_ctz64(word);
            word &= (word - 1u);

            Agent<TSeq> & j = population[j_id];
            if ((j.virus == nullptr) || (j.n_neighbors == 0u))
                continue;

            const size_t s = agent_state[j_id];

            VirusPtr<TSeq> & v = j.virus;

            for (size_t i_id : *j.neighbors)
            {

                // Most neighbors are usually not susceptible; the compact copy
                // of the states rules them out without loading the agent.
                const size_t t = agent_state[i_id];

                if (!push_pushable[t] || push_excluded[t * ns + s])
                    continue;

                Agent<TSeq> & i = population[i_id];

                // An agent with a virus is not susceptible (pulling would
                // refuse it; update_state() reports it).
                if (i.virus != nullptr)
                    continue;

                // Pulling only updates queued agents.
                if (use_queuing && (queue[i_id] <= 0))
                    continue;

                // Exactly the expression the pull uses, in the same order.
                epiworld_double p =
                    (1.0 - i.get_susceptibility_reduction(v, *this)) *
                    v->get_prob_infecting(this) *
                    (1.0 - j.get_transmission_reduction(v, *this))
                    ;

                // No chance of transmission (also catches NaN)
                if (!(p > 0.0))
                    continue;

                int & slot = push_slot[i_id];
                if (slot < 0)
                {

                    slot = static_cast< int >(push_targets.size());
                    push_targets.push_back({i_id, 0.0, 0u, nullptr});

                    #ifdef EPI_DEBUG
                    if (push_default[t])
                        db.n_transmissions_potential++;
                    #endif

                }

                PushTarget & target = push_targets[static_cast< size_t >(slot)];

                // Certain transmission: uniform reservoir among these
                if (p >= 1.0)
                {

                    if (
                        (++target.n_certain == 1u) ||
                        (runif() * static_cast< double >(target.n_certain) < 1.0)
                    )
                        target.candidate = &(*v);

                    continue;

                }

                // A certain transmission wins outright
                if (target.n_certain > 0u)
                    continue;

                // Weighted reservoir: keep this infector with probability
                // r / (sum of odds so far).
                const double odds = static_cast< double >(p) /
                    (1.0 - static_cast< double >(p));
                target.odds += odds;

                if (
                    (target.candidate == nullptr) ||
                    (runif() * target.odds < odds)
                )
                    target.candidate = &(*v);

            }

        }

    }

    // Phase 2: one draw per touched agent. P(no infection) = 1 / (1 + R).
    for (auto & target : push_targets)
    {

        push_slot[target.id] = -1;

        if (target.candidate == nullptr)
            continue;

        if (
            (target.n_certain == 0u) &&
            (runif() < 1.0 / (1.0 + target.odds))
        )
            continue;

        Agent<TSeq> & i = population[target.id];

        #ifdef EPI_DEBUG
        if (push_default[i.state])
            db.n_transmissions_today++;
        #endif

        i.set_virus(*this, *target.candidate);

    }

}

template<typename TSeq>
inline void Model<TSeq>::transmission_update_others()
{

    const size_t ns = static_cast< size_t >(nstates);

    // Pulling refuses an agent in a susceptible state that carries a virus;
    // so does pushing, for the agents a pull would have visited.
    for (size_t s = 0u; s < ns; ++s)
    {

        if (!push_pushable[s] || (state_carriers[s] == 0u))
            continue;

        for (size_t id : state_index_members(s))
        {

            const auto & p = population[id];
            if ((p.virus != nullptr) && (!use_queuing || (queue[id] > 0)))
                throw std::logic_error(
                    std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                    std::string("Agent id ") + std::to_string(p.get_id()) +
                    std::string(" has a virus.")
                    );

        }

    }

    // Everyone in a state with an update function, other than the pushed ones.
    // The index gives them directly, so neither the population nor the queue
    // (which holds every neighbor of every carrier) needs to be scanned. They
    // are marked in a bitset and visited in ascending id order -- the same
    // order whether queuing is on or off -- in O(agents + N / 64), no sort.
    const size_t nwords = (population.size() + 63u) / 64u;
    if (push_visit.size() != nwords)
        push_visit.assign(nwords, 0u);

    for (size_t s = 0u; s < ns; ++s)
        if (state_fun[s] && !push_pushable[s])
            for (size_t id : state_index_members(s))
                push_visit[id >> 6] |= (uint64_t(1) << (id & 63u));

    for (size_t w = 0u; w < nwords; ++w)
    {

        uint64_t word = push_visit[w];
        push_visit[w] = 0u;

        while (word != 0u)
        {

            const size_t id = (w << 6) + epi_ctz64(word);
            word &= (word - 1u);

            // Queued agents only, read as the loop reaches them (a state
            // function may change the queue by editing ties).
            if (use_queuing && (queue[id] <= 0))
                continue;

            auto & p = population[id];
            state_fun[p.state](&p, this);

        }

    }

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::set_transmission_mode(
    TransmissionMode mode,
    double kappa
)
{

    if (!(kappa >= 0.0) || std::isinf(kappa))
        throw std::range_error(
            "The transmission kappa must be a finite, non-negative number."
        );

    transmission_mode = mode;
    transmission_kappa = kappa;
    return *this;

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::set_transmission_mode(
    std::string_view mode,
    double kappa
)
{

    if (mode == "auto")
        return set_transmission_mode(TransmissionMode::automatic, kappa);
    else if (mode == "push")
        return set_transmission_mode(TransmissionMode::push, kappa);
    else if (mode == "pull")
        return set_transmission_mode(TransmissionMode::pull, kappa);

    throw std::invalid_argument(
        "Unknown transmission mode \"" + std::string(mode) +
        "\". Use \"auto\", \"push\", or \"pull\"."
    );

}

template<typename TSeq>
inline TransmissionMode Model<TSeq>::get_transmission_mode() const
{
    return transmission_mode;
}

template<typename TSeq>
inline TransmissionMode Model<TSeq>::get_last_transmission_mode() const
{
    return transmission_mode_last;
}

template<typename TSeq>
inline double Model<TSeq>::get_transmission_kappa() const
{
    return transmission_kappa;
}

#endif
