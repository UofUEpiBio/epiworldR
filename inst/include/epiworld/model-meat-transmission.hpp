#ifndef EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP
#define EPIWORLD_MODEL_MEAT_TRANSMISSION_HPP

/**
 * @file model-meat-transmission.hpp
 * @brief Transmission: from contacts to infections, and network transmission
 * by pushing infection hazards.
 *
 * @details A susceptible agent `i` in a stage with a default infection
 * *pulls*: it scans its contacts, collects the per-contact probabilities
 * `p_ij` of those carrying a virus (in an infectious stage), and one draw
 * (`draw_exposure()`) decides "no infection" or a single infector. The
 * exposures are independent: each has the hazard `h_ij = -log(1 - p_ij)`
 * (times its layer's weight), and with `H_i` the sum of the hazards
 *
 *     P(no infection) = exp(-H_i) = prod_j (1 - p_ij),
 *     P(infected by j) = (h_ij / H_i) * (1 - exp(-H_i)).
 *
 * The draw depends only on the sum of the hazards and on each one's share, so
 * the same outcome can be *pushed*: every agent carrying a virus adds its
 * hazards to its susceptible contacts, each of which keeps a weighted
 * reservoir sample of its infectors; then each touched agent makes one draw.
 * The cost follows the carriers' contacts instead of the susceptibles', which
 * is much cheaper while an outbreak is small.
 *
 * Contacts with `p_ij >= 1` are certain (an infinite hazard): the infection
 * happens and the infector is drawn uniformly among them.
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
    // Decided once per call: several diseases filter the exposures
    return (tracks.size() > 1u) ?
        add_exposures_impl< true >(agent, layer, contacts, n, n_exposures) :
        add_exposures_impl< false >(agent, layer, contacts, n, n_exposures);
}

template< bool Several >
inline size_t Model::add_exposures_impl(
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

        // Several diseases: only those of the draw (see nh_step_infection()),
        // and those the agent can catch
        if (Several && !(nh_trial_mask & v->disease_bit))
            continue;

        const epiworld_double p_inf = nh_infection_prob(
            *agent, *v, Several ? nh_track(v.get()) : tracks[0u]
        );
        if (Several && (p_inf <= 0.0))
            continue;

        const epiworld_double p =
            (1.0 - agent->get_susceptibility_reduction(v, *this, layer_id)) *
            p_inf *
            (1.0 - contact.get_transmission_reduction(v, *this, layer_id))
            ;

        array_double_tmp[n_exposures] = exposure_hazard(p, weight);
        array_virus_tmp[n_exposures] = &(*v);
        array_layer_tmp[n_exposures++] = layer_id;

    }

    return n_exposures;

}

inline epiworld_double Model::exposure_hazard(epiworld_double p, epiworld_double weight)
{
    // A certain exposure (p >= 1) is an infinite hazard; p <= 0 (or NaN) none
    if (p >= 1.0)
        return std::numeric_limits< epiworld_double >::infinity();
    if (!(p > 0.0))
        return 0.0;
    return -weight * std::log1p(-p);
}

inline Virus * Model::draw_exposure(size_t n_exposures)
{

    sampled_layer = -1;

    if (n_exposures == 0u)
        return nullptr;

    // The hazards of the exposures and their sum; certain ones (infinite
    // hazards) are counted in the second half of the scratch array, then one
    // of them is drawn uniformly
    const auto & h = array_double_tmp;
    size_t n_certain = 0u;
    double total = 0.0;
    for (size_t k = 0u; k < n_exposures; ++k)
    {
        if (std::isinf(h[k]))
            ++n_certain;
        else
            total += static_cast< double >(h[k]);
    }

    const double r = runif();
    size_t which = n_exposures;
    if (n_certain > 0u)
    {
        // The (r * n_certain)-th certain exposure
        size_t pick = std::min(
            static_cast< size_t >(r * static_cast< double >(n_certain)),
            n_certain - 1u
        );
        for (which = 0u; which < n_exposures; ++which)
            if (std::isinf(h[which]) && (pick-- == 0u))
                break;
    }
    else
    {
        // No infection with probability exp(-total); otherwise the infector
        // is drawn with probability proportional to its hazard, reusing r
        // (uniform on [p_none, 1) given an infection)
        const double p_none = std::exp(-total);
        if ((total <= 0.0) || (r < p_none))
            return nullptr;

        double target = (r - p_none) / (1.0 - p_none) * total;
        for (which = 0u; which < (n_exposures - 1u); ++which)
        {
            target -= static_cast< double >(h[which]);
            if ((target < 0.0) && (h[which] > 0.0))
                break;
        }

        // Rounding can leave the walk at the last entry; take the last one
        // with a hazard
        while (!(h[which] > 0.0))
            --which;
    }

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
    const epiworld_fast_uint state = agent_state[id];

    size_t n_contacts = 0u;
    size_t n_exposures = 0u;
    for (const auto & layer : layers)
    {

        // A closed layer, or one whose statuses leave the agent out, has no
        // contacts for it (the states out of every layer are the layers' to
        // handle, as before)
        const std::vector< char > * out = layer->own_exclusions();
        if (layer->is_closed() || ((out != nullptr) && (*out)[state]))
            continue;

        layer_contacts.clear();
        const size_t n = layer->sample(id, *this, layer->silent(), layer_contacts);
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

    // A natural history says which stages are pushed to (several diseases
    // are pulled: several-diseases.md, Section 4.2)
    if (tracks.size() == 1u)
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
                tracks[0u].nh->stages[static_cast< size_t >(nh_stage_of(cu))].pushable &&
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
    //
    // A mixing layer costs each visited agent one binomial draw per group it
    // draws from (all the groups with a contact matrix, its own within
    // groups), EPI_TRANSMISSION_DRAW_COST ties each, on either side: pushing
    // draws whom a carrier meets, pulling whom a susceptible meets.
    const TransmissionSums sums = get_transmission_sums();

    double draws = 0.0;
    for (const auto & layer : layers)
        if (const auto * mix = layer->as_mixing())
            draws += mix->within_groups() ?
                1.0 : static_cast< double >(mix->get_n_groups());

    const double per_agent = EPI_TRANSMISSION_AGENT_COST +
        EPI_TRANSMISSION_DRAW_COST * draws;

    const double cost_push =
        static_cast< double >(sums.carrier_degree) +
        per_agent * static_cast< double >(sums.carriers);

    const double cost_pull =
        static_cast< double >(sums.susceptible_degree) +
        per_agent * static_cast< double >(sums.susceptibles);

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
    for (const auto & layer : layers)
        record = record || layer->is_tracing();

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

    // A pull visits only queued agents, unless a layer is blind to the queue
    // (a contact matrix): then it visits everyone, and so does the push
    const bool queue_filter = use_queuing && !has_queue_blind_layer();

    // Phase 1: every carrier that can transmit adds its hazards to the eligible
    // agents it reaches. Nothing changes state until events_run(), so this sees
    // the model as it was at the start of the step, as pulling does.
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

            // Every layer, in the order the layers were added, as a pull reads
            // them
            for (const auto & layer_ptr : layers)
            {

                ContactLayer * layer = layer_ptr.get();

                // A closed layer, or one whose statuses leave the carrier out,
                // reaches no one (a carrier out of every layer is no source)
                const std::vector< char > * out = layer->own_exclusions();
                if (layer->is_closed() || ((out != nullptr) && (*out)[agent_state[j_id]]))
                    continue;

                // Layer rules, as add_exposures() applies them in a pull
                const bool weighted = layer->has_weight();
                const epiworld_double weight = weighted ?
                    layer->get_weight(*this) : 1.0;
                const bool admitted = layer->admits(*v);
                const int layer_id = static_cast< int >(layer->get_id());

                // The agents the carrier reaches through the layer. A
                // network: its ties. Another layer: the agents it meets this
                // step, drawn as likely as their pulls would draw it (see
                // ContactLayer::sample_reverse()). One loop over both keeps
                // the body inlined.
                const size_t * reached;
                size_t n_reached;
                if (auto * net = layer->as_network())
                {
                    reached = net->ties[j_id].data();
                    n_reached = net->ties[j_id].size();
                }
                else
                {
                    push_reached.clear();
                    layer->sample_reverse(j_id, *this, push_reached);
                    reached = push_reached.data();
                    n_reached = push_reached.size();
                }

                for (size_t r = 0u; r < n_reached; ++r)
                {

                    const size_t i_id = reached[r];

                    // Most agents reached are usually not susceptible; the
                    // compact copy of the states rules them out without
                    // loading the agent.
                    const size_t t = agent_state[i_id];

                    if (!push_pushable[t] || ((out != nullptr) && (*out)[t]))
                        continue;

                    Agent & i = population[i_id];

                    // An agent with a virus is not susceptible (pulling would
                    // refuse it; update_state() reports it).
                    if (i.virus != nullptr)
                        continue;

                    // Pulling only updates queued agents.
                    if (queue_filter && (queue[i_id] <= 0))
                        continue;

                    // The post-sampling callback and the layer's log see every
                    // eligible contact, also the ones that cannot transmit, as
                    // a pull reports them. The carriers are visited in
                    // ascending id order, so the pairs come out grouped by
                    // infectious agent.
                    if constexpr (Record)
                    {
                        if (post_sampling_on)
                            register_sampled_contact(j_id, i_id);
                        layer->trace(j_id, i_id);
                    }

                    // A layer of weight 0, or closed to this virus, transmits
                    // nothing (its contacts were still reported above)
                    if (!admitted || (weighted && (weight <= 0.0)))
                        continue;

                    // Exactly the expression the pull uses, in the same order.
                    const epiworld_double h = exposure_hazard(
                        (1.0 - i.get_susceptibility_reduction(v, *this, layer_id)) *
                        nh_infection_prob(i, *v, tracks[0u]) *
                        (1.0 - j.get_transmission_reduction(v, *this, layer_id)),
                        weight
                    );

                    // No chance of transmission (also catches NaN)
                    if (!(h > 0.0))
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
                    if (std::isinf(h))
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
                    // h / (sum of hazards so far).
                    target.hazard += static_cast< double >(h);

                    if (
                        (target.candidate == nullptr) ||
                        (runif() * target.hazard < static_cast< double >(h))
                    )
                    {
                        target.candidate = &(*v);
                        target.layer = layer_id;
                    }

                }

            }

        }

    }

    // Phase 2: one draw per touched agent. P(no infection) = exp(-H).
    for (auto & target : push_targets)
    {

        push_slot[target.id] = -1;

        if (target.candidate == nullptr)
            continue;

        if (
            (target.n_certain == 0u) &&
            (runif() < std::exp(-target.hazard))
        )
            continue;

        Agent & i = population[target.id];

        #ifdef EPI_DEBUG
        db.n_transmissions_today++;
        #endif

        sampled_layer = target.layer;
        if (!tracks.empty())
            nh_propose(
                tracks[0u],
                tracks[0u].nh->stages[static_cast< size_t >(nh_stage_of(i.state))].infection,
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
            if ((p.virus != nullptr) && (!use_queuing || has_queue_blind_layer() || (queue[id] > 0)))
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
            // function may change the queue by editing ties), unless a layer
            // is blind to the queue: a pull then visits everyone
            if (use_queuing && !has_queue_blind_layer() && (queue[id] <= 0))
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
