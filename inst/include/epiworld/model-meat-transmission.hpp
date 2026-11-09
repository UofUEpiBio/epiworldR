#ifndef EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP
#define EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP

/**
 * @file model-meat-transmission.hpp
 * @brief Transmission: from contacts to infections, and network transmission
 * by pushing infection odds.
 *
 * @details A susceptible agent `i` in a stage with a default infection
 * *pulls*: it scans its neighbors, collects the per-contact probabilities
 * `p_ij` of those carrying a virus (in an infectious stage), and `roulette()` draws "no infection" or a
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

inline void Model::state_masks_build()
{

    const size_t ns = static_cast< size_t >(nstates);

    auto check = [ns](int s, const char * what) {
        if ((s < 0) || (static_cast< size_t >(s) >= ns))
            throw std::range_error(
                std::string("The ") + what + " state " + std::to_string(s) +
                " is out of range (the model has " + std::to_string(ns) +
                " states)."
            );
    };

    state_excluded.assign(ns, 0);
    for (int s : excluded_states)
    {
        check(s, "excluded");
        state_excluded[s] = 1;
    }

    // An excluded agent does not transmit either
    state_not_infectious.clear();
    if (infectious_states.empty() && excluded_states.empty())
        return;

    state_not_infectious.assign(ns, !infectious_states.empty());
    for (int s : infectious_states)
    {
        check(s, "infectious");
        state_not_infectious[s] = false;
    }

    for (int s : excluded_states)
        state_not_infectious[s] = true;

}

inline size_t Model::add_exposures(
    Agent * agent,
    const ContactLayer * layer,
    const size_t * contacts,
    size_t n,
    size_t n_exposures
)
{

    const bool weighted = (layer != nullptr) && layer->has_weight();
    const epiworld_double weight = weighted ? layer->get_weight(*this) : 1.0;
    const int layer_id = (layer == nullptr) ? -1 : static_cast< int >(layer->get_id());

    // A layer of weight 0 transmits nothing
    if (weighted && (weight <= 0.0))
        return n_exposures;

    reserve_tmp(n_exposures + n);

    for (size_t k = 0u; k < n; ++k)
    {

        auto & contact = population[contacts[k]];
        auto & v = contact.get_virus();
        if (v == nullptr)
            continue;

        if ((layer != nullptr) && !layer->admits(*v))
            continue;

        epiworld_double p =
            (1.0 - agent->get_susceptibility_reduction(v, *this)) *
            nh_infection_prob(*agent, *v) *
            (1.0 - contact.get_transmission_reduction(v, *this))
            ;

        // The weight scales the odds p / (1 - p); a certain exposure stays so
        if (weighted && (p < 1.0))
            p = weight * p / (1.0 - p + weight * p);

        array_double_tmp[n_exposures] = p;
        array_virus_tmp[n_exposures] = &(*v);
        array_layer_tmp[n_exposures++] = layer_id;

    }

    return n_exposures;

}

inline Virus * Model::draw_exposure(size_t n_exposures)
{

    sampled_layer = -1;

    if (n_exposures == 0u)
        return nullptr;

    const int which = roulette(n_exposures, this);

    if (which < 0)
        return nullptr;

    sampled_layer = array_layer_tmp[which];
    return array_virus_tmp[which];

}

inline Virus * Model::sample_virus(
    Agent * agent,
    const size_t * contacts,
    size_t n,
    const ContactLayer * layer
)
{

    return draw_exposure(add_exposures(agent, layer, contacts, n, 0u));

}

inline Virus * Model::sample_virus(Agent * agent)
{

    if (agent->get_virus() != nullptr)
        throw std::logic_error(
            "Sampling a virus for agent " + std::to_string(agent->get_id()) +
            ", which already has one, makes no sense."
        );

    const size_t id = static_cast< size_t >(agent->get_id());
    const std::vector< bool > * exclude =
        state_not_infectious.empty() ? nullptr : &state_not_infectious;

    size_t n_contacts = 0u;
    size_t n_exposures = 0u;
    for (const auto & layer : layers)
    {

        layer_contacts.clear();
        const size_t n = layer->sample(id, *this, exclude, layer_contacts);
        if (n == 0u)
            continue;

        n_contacts += n;
        n_exposures = add_exposures(
            agent, layer.get(), layer_contacts.data(), n, n_exposures
        );

    }

    if (n_contacts == 0u)
    {
        sampled_layer = -1;
        return nullptr;
    }

    #ifdef EPI_DEBUG
    db.n_transmissions_potential++;
    #endif

    Virus * virus = draw_exposure(n_exposures);

    #ifdef EPI_DEBUG
    if (virus != nullptr)
        db.n_transmissions_today++;
    #endif

    return virus;

}

inline bool Model::transmission_prepare()
{

    const size_t ns = static_cast< size_t >(nstates);

    push_pushable.assign(ns, 0);
    push_source_ok.assign(ns, 0);

    bool any = false;

    // A natural history says which stages are pushed to
    if (nh)
    {
        // A cell is pushed to if its stage's only drawn exit is a default
        // infection, its status has no drawn exit, and it is in the layers
        for (size_t c = 0u; c < ns; ++c)
        {
            const unsigned int cu = static_cast< unsigned int >(c);
            const bool status_drawn = st &&
                st->stages[static_cast< size_t >(nh_status_of(cu))].needs_visit;
            const bool excluded = !state_excluded.empty() && state_excluded[c];
            if (
                nh->stages[static_cast< size_t >(nh_stage_of(cu))].pushable &&
                !status_drawn && !excluded
            )
            {
                push_pushable[c] = 1;
                any = true;
            }
        }
    }

    // The infectious states are the sources, if anyone is pushed to
    if (any)
        for (size_t s = 0u; s < ns; ++s)
            push_source_ok[s] =
                state_not_infectious.empty() || !state_not_infectious[s];

    return any;

}

inline bool Model::transmission_choose_push() const
{

    if (transmission_mode == TransmissionMode::pull)
        return false;

    if (transmission_mode == TransmissionMode::push)
        return true;

    // Pushing visits every carrier that can transmit and walks its ties;
    // pulling visits every susceptible agent and walks its ties. Each visit
    // also costs a fixed amount, EPI_TRANSMISSION_AGENT_COST ties: at the low
    // degrees of contact networks that is comparable to the ties themselves,
    // so it matters which side has more agents. The sums are kept per state, so
    // this is O(number of states). It deliberately ignores the queue: the
    // decision -- and so the random stream -- is the same with queuing on or
    // off. kappa < 1 makes pushing harder to choose: the queue spares a pull
    // the susceptibles with no infectious neighbor, which the sums do not see,
    // and the per-agent term raises the pull side wherever susceptibles
    // outnumber carriers. The default was fitted so the two balance.
    const TransmissionSums sums = get_transmission_sums();

    const double cost_push =
        static_cast< double >(sums.carrier_degree) +
        EPI_TRANSMISSION_AGENT_COST * static_cast< double >(sums.carriers);

    const double cost_pull =
        static_cast< double >(sums.susceptible_degree) +
        EPI_TRANSMISSION_AGENT_COST * static_cast< double >(sums.susceptibles);

    return cost_push <= transmission_kappa * cost_pull;

}

inline typename Model::TransmissionSums
Model::get_transmission_sums() const
{

    TransmissionSums sums = {0u, 0u, 0u, 0u};

    if (!state_index_ready || (push_pushable.size() != static_cast< size_t >(nstates)))
        return sums;

    for (size_t s = 0u; s < static_cast< size_t >(nstates); ++s)
    {
        if (push_source_ok[s])
        {
            sums.carrier_degree += state_carrier_degree[s];
            sums.carriers += state_carriers[s];
        }
        if (push_pushable[s])
        {
            sums.susceptible_degree += state_degree[s];
            sums.susceptibles += state_start[s + 1u] - state_start[s];
        }
    }

    return sums;

}

inline void Model::transmission_push()
{

    // The recording path is chosen once per step, not per contact
    bool record = post_sampling_on;
    for (const auto * net : network_layers)
        record = record || net->is_tracing();

    if (record)
        transmission_push_impl<true>();
    else
        transmission_push_impl<false>();

}

template<bool Record>
inline void Model::transmission_push_impl()
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

            Agent & j = population[j_id];
            if (j.virus == nullptr)
                continue;

            VirusPtr & v = j.virus;

            // Every network, in the order the layers were added, as a pull
            // reads them
            for (auto * net : network_layers)
            {

                // Layer rules, as add_exposures() applies them in a pull
                const bool weighted = net->has_weight();
                const epiworld_double weight = weighted ?
                    net->get_weight(*this) : 1.0;
                const bool admitted = net->admits(*v);
                const int layer_id = static_cast< int >(net->get_id());

                for (size_t i_id : net->ties[j_id])
                {

                    // Most neighbors are usually not susceptible; the compact copy
                    // of the states rules them out without loading the agent.
                    const size_t t = agent_state[i_id];

                    if (!push_pushable[t])
                        continue;

                    Agent & i = population[i_id];

                    // An agent with a virus is not susceptible (pulling would
                    // refuse it; update_state() reports it).
                    if (i.virus != nullptr)
                        continue;

                    // Pulling only updates queued agents.
                    if (use_queuing && (queue[i_id] <= 0))
                        continue;

                    // The post-sampling callback and the layer's log see every
                    // eligible contact, also the ones that cannot transmit, as a
                    // pull reports them. The carriers are visited in ascending id
                    // order, so the pairs come out grouped by infectious agent.
                    if constexpr (Record)
                    {
                        if (post_sampling_on)
                            register_sampled_contact(j_id, i_id);
                        net->trace(j_id, i_id);
                    }

                    // A layer of weight 0, or closed to this virus, transmits
                    // nothing (its contacts were still reported above)
                    if (!admitted || (weighted && (weight <= 0.0)))
                        continue;

                    // Exactly the expression the pull uses, in the same order.
                    epiworld_double p =
                        (1.0 - i.get_susceptibility_reduction(v, *this)) *
                        nh_infection_prob(i, *v) *
                        (1.0 - j.get_transmission_reduction(v, *this))
                        ;

                    // The weight scales the odds; a certain exposure stays so
                    if (weighted && (p < 1.0))
                        p = weight * p / (1.0 - p + weight * p);

                    // No chance of transmission (also catches NaN)
                    if (!(p > 0.0))
                        continue;

                    int & slot = push_slot[i_id];
                    if (slot < 0)
                    {

                        slot = static_cast< int >(push_targets.size());
                        push_targets.push_back({i_id, 0.0, 0u, nullptr, -1});

                        #ifdef EPI_DEBUG
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
                        {
                            target.candidate = &(*v);
                            target.layer = layer_id;
                        }

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
                    {
                        target.candidate = &(*v);
                        target.layer = layer_id;
                    }

                }

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

        Agent & i = population[target.id];

        #ifdef EPI_DEBUG
        db.n_transmissions_today++;
        #endif

        sampled_layer = target.layer;
        if (nh)
            nh_propose(
                nh->stages[static_cast< size_t >(nh_stage_of(i.state))].infection,
                i, target.candidate
            );
        else
            i.set_virus(*this, *target.candidate);

    }

}

inline void Model::transmission_update_others()
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
                    std::string("A pushed-to stage has agents WITH viruses! ") +
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

inline Model & Model::set_transmission_mode(
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

inline Model & Model::set_transmission_mode(
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

inline TransmissionMode Model::get_transmission_mode() const
{
    return transmission_mode;
}

inline TransmissionMode Model::get_last_transmission_mode() const
{
    return transmission_mode_last;
}

inline double Model::get_transmission_kappa() const
{
    return transmission_kappa;
}

#endif
