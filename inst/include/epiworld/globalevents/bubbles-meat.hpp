#ifndef EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP
#define EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP

// Standard library headers are included at global scope by epiworld.hpp (this
// file is included from within `namespace epiworld`, so system headers must not
// be re-included here).
#include "bubbles-bones.hpp"

template<typename TSeq>
inline BubbleTool<TSeq>::BubbleTool(
    std::string name,
    std::string event_name
) : Tool<TSeq>(name), _event_name(std::move(event_name))
{
}

template<typename TSeq>
inline epiworld_double BubbleTool<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> * model
)
{

    // Binding to the policy of *this* model, once. clone_ptr() clears the
    // pointer, so a copy of this tool -- in another agent, or in a copy of the
    // model -- resolves against its own model rather than inheriting ours.
    if (_policy == nullptr)
    {

        _policy = Bubbles<TSeq>::get_from(*model, _event_name);

        if (_policy == nullptr)
            throw std::logic_error(
                "BubbleTool: the intervention '" + _event_name +
                "' is not installed on this model."
            );

    }

    return _policy->susceptibility_reduction(
        this->get_agent(), v->get_agent(), model
    );

}

template<typename TSeq>
inline std::unique_ptr<Tool<TSeq>> BubbleTool<TSeq>::clone_ptr() const
{
    auto ans = std::make_unique<BubbleTool<TSeq>>(*this);
    ans->_policy = nullptr; // the copy resolves against its own model
    return ans;
}

template<typename TSeq>
inline Bubbles<TSeq>::Bubbles(
    std::vector< size_t > household_id,
    BubbleFlavor flavor,
    size_t group_size,
    epiworld_double transmission_factor,
    int start_day,
    int end_day,
    int rewire_every,
    std::string name,
    size_t max_households,
    std::string param_name,
    BubbleTies ties
) :
    household_id(std::move(household_id)),
    flavor(flavor),
    group_size(group_size),
    max_households(max_households),
    transmission_factor(transmission_factor),
    start_day(start_day),
    end_day(end_day),
    rewire_every(rewire_every),
    param_name(std::move(param_name)),
    ties(ties)
{

    this->set_name(name);
    this->set_day(-99); // runs at the end of every day

    if ((this->transmission_factor < 0.0) || (this->transmission_factor > 1.0))
        throw std::range_error(
            "Bubbles: transmission_factor must be in [0, 1]."
        );

    if ((flavor == BubbleFlavor::Household) && (group_size < 1u))
        throw std::range_error(
            "Bubbles: group_size (households per bubble) must be >= 1."
        );

    if ((flavor == BubbleFlavor::Peer) && (max_households < 2u))
        throw std::range_error(
            "Bubbles: max_households must be >= 2 for the Peer flavor."
        );

    if ((this->end_day >= 0) && (this->end_day <= this->start_day))
        throw std::range_error(
            "Bubbles: end_day must be greater than start_day (or negative)."
        );

}

template<typename TSeq>
inline Bubbles<TSeq>::HeldTies::HeldTies(
    Model<TSeq> * model,
    const Bubbles<TSeq> & self
)
{

    // Every policy's books, this one's included.
    std::vector< const Bubbles<TSeq> * > holders = {&self};
    for (size_t e = 0u; e < model->get_n_globalevents(); ++e)
    {

        const auto * other = dynamic_cast< const Bubbles<TSeq> * >(
            &model->get_globalevent(e)
        );

        if ((other != nullptr) && (other != &self))
            holders.push_back(other);

    }

    // A policy that has not been set up for this model yet may still hold
    // ties from another network; only ties between this model's agents count.
    size_t n = model->size();
    std::vector< std::pair< size_t, size_t > > ties;
    for (const auto * holder : holders)
        for (const auto & tie : holder->created_ties)
            if ((tie.first < n) && (tie.second < n))
                ties.push_back(tie);

    if (ties.empty())
        return;

    // Compressed rows: the held partners of agent `a` are
    // partner[start[a]] .. partner[start[a + 1] - 1].
    start.assign(n + 1u, 0u);
    for (const auto & tie : ties)
    {
        ++start[tie.first + 1u];
        ++start[tie.second + 1u];
    }

    for (size_t a = 0u; a < n; ++a)
        start[a + 1u] += start[a];

    partner.resize(start[n]);
    std::vector< size_t > fill(start.begin(), start.end() - 1);
    for (const auto & tie : ties)
    {
        partner[fill[tie.first]++]  = tie.second;
        partner[fill[tie.second]++] = tie.first;
    }

    marked.assign(n, 0);

}

template<typename TSeq>
inline void Bubbles<TSeq>::HeldTies::focus(size_t a)
{

    if (marked.empty())
        return;

    for (size_t k = start[focused]; k < start[focused + 1u]; ++k)
        marked[partner[k]] = 0;

    focused = a;

    for (size_t k = start[a]; k < start[a + 1u]; ++k)
        marked[partner[k]] = 1;

}

template<typename TSeq>
inline bool Bubbles<TSeq>::HeldTies::contains(size_t b) const
{
    return !marked.empty() && (marked[b] != 0);
}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_household(
    Model<TSeq> * model,
    HeldTies & held
)
{

    // Map household label -> compact index, and list unique households.
    std::unordered_map< size_t, size_t > hh_index;
    std::vector< size_t > hh_labels;
    for (size_t a = 0u; a < household_id.size(); ++a)
    {
        size_t h = household_id[a];
        if (hh_index.find(h) == hh_index.end())
        {
            hh_index[h] = hh_labels.size();
            hh_labels.push_back(h);
        }
    }

    size_t nh = hh_labels.size();

    // Build the household contact graph: households h1 and h2 are adjacent when
    // at least one member of h1 is connected to a member of h2 in the contact
    // network. Bubbles are grown along these ties.
    //
    // Grouping households that share NO tie would be a no-op: the intervention
    // can only suppress transmission along existing edges, never create new
    // ones, so bubbling two unconnected households changes nothing. Pairing at
    // random therefore degenerates to the household-only lockdown. It also
    // matches the policy being modelled: a household picks a bubble partner it
    // actually socialises with.
    //
    // Ties a bubble policy is holding are not contacts anybody has, so they do
    // not count (see compute_partition()).
    std::vector< std::vector< size_t > > hh_adj(nh);
    auto & pop_all = model->get_agents();
    for (size_t a = 0u; a < household_id.size(); ++a)
    {
        size_t ha = hh_index[household_id[a]];
        held.focus(a);
        for (auto * nb : pop_all[a].get_neighbors(*model))
        {
            size_t b = static_cast< size_t >(nb->get_id());
            if (household_id[a] == household_id[b])
                continue;
            if (held.contains(b))
                continue;
            hh_adj[ha].push_back(hh_index[household_id[b]]);
        }
    }

    // De-duplicate each adjacency list.
    for (auto & adj : hh_adj)
    {
        std::sort(adj.begin(), adj.end());
        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
    }

    // Shuffle household order (Fisher-Yates with the model RNG).
    std::vector< size_t > order(nh);
    for (size_t i = 0u; i < nh; ++i)
        order[i] = i;

    for (size_t i = nh; i > 1u; --i)
    {
        size_t j = static_cast< size_t >(model->runif_index(static_cast<uint32_t>(i)));
        std::swap(order[i - 1u], order[j]);
    }

    // Grow each bubble from a seed household by repeatedly absorbing a random
    // household that is *connected* to the bubble, up to `group_size`
    // households. A household with no unassigned connected candidates simply
    // ends up in a smaller bubble (possibly alone) -- you can only bubble with
    // someone you already have contact with.
    std::vector< int > hh_bubble(nh, -1);
    std::vector< size_t > candidates;
    int next_bubble = 0;

    for (size_t pos = 0u; pos < nh; ++pos)
    {

        size_t seed = order[pos];
        if (hh_bubble[seed] != -1)
            continue;

        int b = next_bubble++;
        hh_bubble[seed] = b;
        size_t members = 1u;

        // Frontier of households connected to the bubble. May contain stale
        // (already assigned) or repeated entries; repeats make a household that
        // is tied to several members proportionally more likely to be picked.
        candidates.clear();
        for (size_t x : hh_adj[seed])
            if (hh_bubble[x] == -1)
                candidates.push_back(x);

        while ((members < group_size) && !candidates.empty())
        {

            size_t idx = static_cast< size_t >(
                model->runif_index(static_cast<uint32_t>(candidates.size()))
            );
            size_t pick = candidates[idx];
            candidates[idx] = candidates.back();
            candidates.pop_back();

            if (hh_bubble[pick] != -1) // stale entry
                continue;

            hh_bubble[pick] = b;
            ++members;

            for (size_t x : hh_adj[pick])
                if (hh_bubble[x] == -1)
                    candidates.push_back(x);

        }

    }

    // Assign each agent the bubble of its household.
    for (size_t a = 0u; a < household_id.size(); ++a)
        bubble_id[a] = hh_bubble[hh_index[household_id[a]]];

}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_peer(
    Model<TSeq> * model,
    HeldTies & held
)
{

    size_t n = household_id.size();

    // Map household label -> compact index.
    std::unordered_map< size_t, size_t > hh_index;
    std::vector< size_t > hh_labels;
    for (size_t a = 0u; a < n; ++a)
    {
        size_t h = household_id[a];
        if (hh_index.find(h) == hh_index.end())
        {
            hh_index[h] = hh_labels.size();
            hh_labels.push_back(h);
        }
    }

    size_t nh = hh_labels.size();

    // Disjoint-set (union-find) over households, tracking the number of
    // households in each set so bubbles can be capped.
    std::vector< size_t > parent(nh), set_size(nh, 1u);
    for (size_t i = 0u; i < nh; ++i)
        parent[i] = i;

    auto find = [&parent](size_t x) -> size_t {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]]; // path halving
            x = parent[x];
        }
        return x;
    };

    // Agents choose in random order, each drawing peers from the contacts that
    // are still available. A household whose bubble is full drops out of the
    // pool: any choice involving it is declined, and its own members stop
    // choosing. This cap is what makes the policy's exclusivity bite -- without
    // it the merges percolate, the household graph becomes connected, and every
    // household ends up in one giant bubble, imposing no restriction at all.
    std::vector< size_t > agent_order(n);
    for (size_t i = 0u; i < n; ++i)
        agent_order[i] = i;

    for (size_t i = n; i > 1u; --i)
    {
        size_t j = static_cast< size_t >(
            model->runif_index(static_cast<uint32_t>(i))
        );
        std::swap(agent_order[i - 1u], agent_order[j]);
    }

    auto & pop = model->get_agents();
    std::vector< size_t > ext;

    for (size_t oi = 0u; oi < n; ++oi)
    {

        size_t a  = agent_order[oi];
        size_t ha = hh_index[household_id[a]];

        // This agent's household is already in a full bubble: it is out of the
        // pool and cannot take anyone else in.
        if (set_size[find(ha)] >= max_households)
            continue;

        // Households of this agent's contacts outside its own household. Ties
        // a bubble policy is holding are not contacts (see
        // compute_partition()).
        ext.clear();
        held.focus(a);
        for (auto * nb : pop[a].get_neighbors(*model))
        {
            size_t nid = static_cast< size_t >(nb->get_id());
            if ((household_id[nid] != household_id[a]) && !held.contains(nid))
                ext.push_back(hh_index[household_id[nid]]);
        }

        // Keep drawing until the agent has made `group_size` choices or no
        // contact is left that its bubble can still take in.
        size_t chosen = 0u;
        while ((chosen < group_size) && !ext.empty())
        {

            size_t idx = static_cast< size_t >(
                model->runif_index(static_cast<uint32_t>(ext.size()))
            );
            size_t hb = ext[idx];
            ext[idx] = ext.back();
            ext.pop_back();

            size_t ra = find(ha);
            size_t rb = find(hb);

            if (ra == rb) // already sharing a bubble
                continue;

            if ((set_size[ra] + set_size[rb]) > max_households)
                continue; // that bubble is full: not available

            parent[ra] = rb;
            set_size[rb] += set_size[ra];
            ++chosen;

            if (set_size[rb] >= max_households)
                break; // this bubble is now full

        }

    }

    // Compact the component roots to 0..K-1 and label agents.
    std::unordered_map< size_t, int > root_label;
    int next_label = 0;
    for (size_t a = 0u; a < n; ++a)
    {
        size_t root = find(hh_index[household_id[a]]);
        auto it = root_label.find(root);
        if (it == root_label.end())
        {
            root_label[root] = next_label;
            bubble_id[a] = next_label;
            ++next_label;
        }
        else
        {
            bubble_id[a] = it->second;
        }
    }

}

template<typename TSeq>
inline bool Bubbles<TSeq>::wants_ties(Model<TSeq> * model) const
{

    if (ties != BubbleTies::Complete)
        return false;

    // Global events run *after* the day's transitions, so the network this
    // leaves behind is the one the next step will use -- which is also when the
    // tool starts (or stops) damping. Asking about `today() + 1` keeps the ties
    // and the damping switching on the same day.
    int next = model->today() + 1;

    // ... and when there is no next step, nothing needs ties. Withdrawing them
    // on the last day is what keeps the intervention from outliving its run:
    // Model::run() takes no population backup, so ties left in the network
    // would still be there when the model is run again, or would be captured by
    // the backup run_multiple() takes. `ndays == 0` means the day loop is being
    // driven by hand and there is no known end, so the question does not apply.
    size_t ndays = static_cast< size_t >(model->get_ndays());
    if ((ndays > 0u) && (next > static_cast< int >(ndays)))
        return false;

    return is_active(next);

}

template<typename TSeq>
inline void Bubbles<TSeq>::build_ties(Model<TSeq> * model)
{

    // Agents grouped by bubble. The labels compute_partition() hands out are
    // consecutive from zero, so this is an indexed bucket rather than a hash
    // map: the order ties are created in decides the order they take in each
    // agent's neighbor list, and that decides which transmitter roulette()
    // picks. A container with unspecified iteration order would make runs
    // depend on the standard library rather than on the seed.
    int n_bubbles = 0;
    for (int b : bubble_id)
        if (b >= n_bubbles)
            n_bubbles = b + 1;

    std::vector< std::vector< size_t > > members(
        static_cast< size_t >(n_bubbles)
    );

    for (size_t a = 0u; a < bubble_id.size(); ++a)
        if (bubble_id[a] >= 0)
            members[static_cast< size_t >(bubble_id[a])].push_back(a);

    // What the sampler can take is a *degree*, not a bubble size: when pulling,
    // roulette() uses two slots per candidate in the fixed scratch array, so an
    // agent may have at most `array_double_tmp.size() / 2` infectious neighbors
    // before it throws. Pushing has no such ceiling, but the automatic mode may
    // pull on any day, so the check does not depend on the mode. Completing a
    // bubble adds ties on top of the ones an agent already has outside it, so
    // the number that matters is what each member's degree will be afterwards
    // -- a bubble small enough to look harmless can still push a
    // well-connected member over.
    //
    // The policy answers for the ties it adds, not for the network it was
    // handed. A member the network already put past the ceiling -- a hub -- is
    // fine as long as its bubble adds nothing to it; what is refused is a
    // clique that takes a member past the ceiling, or adds to one already
    // there.
    //
    // Checked for every member before a single tie is created, so a bubble that
    // is too large fails without leaving the network half-rewritten.
    size_t max_degree = model->array_double_tmp.size() / 2u;

    for (auto & who : members)
    {

        if (who.size() < 2u)
            continue;

        for (size_t x : who)
        {

            // Ties to bubble-mates the agent already has are not added twice.
            size_t already = 0u;
            for (auto * nb : model->get_agent(x).neighbors_view(*model))
                if (bubble_id[static_cast< size_t >(nb->get_id())] ==
                    bubble_id[x])
                    ++already;

            size_t degree    = model->get_agent(x).get_n_neighbors();
            size_t projected = degree + (who.size() - 1u) - already;

            if ((projected > max_degree) && (projected > degree))
                throw std::length_error(
                    "Bubbles: completing a bubble of " +
                    std::to_string(who.size()) + " would give agent " +
                    std::to_string(x) + " a degree of " +
                    std::to_string(projected) + " (up from " +
                    std::to_string(degree) + "), above the " +
                    std::to_string(max_degree) + " neighbors the virus sampler "
                    "can weigh. Reduce group_size or max_households."
                );

        }

    }

    for (auto & who : members)
    {

        for (size_t i = 0u; i < who.size(); ++i)
        {
            for (size_t j = i + 1u; j < who.size(); ++j)
            {

                // Only ties this actually created are recorded: the ones the
                // network already had are not the intervention's to withdraw.
                if (model->add_edge(who[i], who[j]))
                    created_ties.emplace_back(who[i], who[j]);

            }
        }

    }

    ties_epoch = last_epoch;

}

template<typename TSeq>
inline Bubbles<TSeq> * Bubbles<TSeq>::heir_of(
    Model<TSeq> * model,
    size_t i,
    size_t j
)
{

    // Nothing to hand over unless somebody else is holding a bubble open.
    for (size_t e = 0u; e < model->get_n_globalevents(); ++e)
    {

        auto * other =
            dynamic_cast< Bubbles<TSeq> * >(&model->get_globalevent(e));

        if ((other == nullptr) || (other == this))
            continue;

        // Only a Complete policy whose bubble is about to be in force.
        if (!other->wants_ties(model))
            continue;

        const auto & other_id = other->bubble_id;
        if ((i >= other_id.size()) || (j >= other_id.size()))
            continue;

        if ((other_id[i] < 0) || (other_id[i] != other_id[j]))
            continue;

        return other;

    }

    return nullptr;

}

template<typename TSeq>
inline void Bubbles<TSeq>::withdraw_ties(
    Model<TSeq> * model,
    bool hand_over
)
{

    for (auto & tie : created_ties)
    {

        // Two policies can want the same tie -- household bubbles and school
        // bubbles, say -- but only the one that happened to create it has it on
        // its books. Dropping it here would take it away from a bubble that is
        // still open, so it is handed to that policy instead of being removed.
        //
        // That holds even for a tie something else has taken away in the
        // meantime -- an event isolating an agent, say. The tie is still the
        // bubbles', only suspended: if the isolation ends while the heir's
        // bubble is up, the tie comes back, and it is the heir's to withdraw.
        if (hand_over)
        {

            Bubbles<TSeq> * heir = heir_of(model, tie.first, tie.second);

            if (heir != nullptr)
            {
                heir->created_ties.push_back(tie);
                continue;
            }

        }

        // A tie that is not there -- taken away, and never put back -- is
        // simply skipped: rm_edge() does nothing, and it leaves the books.
        model->rm_edge(tie.first, tie.second);

    }

    created_ties.clear();
    ties_epoch = -1;

}

template<typename TSeq>
inline void Bubbles<TSeq>::restore_network(Model<TSeq> * model)
{
    withdraw_ties(model, false);
}

template<typename TSeq>
inline void Bubbles<TSeq>::sync_ties(Model<TSeq> * model)
{

    bool up = (ties_epoch >= 0);

    if (!wants_ties(model))
    {

        if (up)
            withdraw_ties(model, true);

        return;

    }

    // Standing: the ordinary day, and nothing to do. A standing clique always
    // belongs to the partition in force, since the daily event withdraws it
    // before drawing a new one. It is not re-checked either. A tie that something else took away
    // -- an event isolating an agent, say -- was taken deliberately, and
    // putting it back would undo that. It stays on the books, so that if it is
    // restored while the bubble is up it still comes down with the bubble. The
    // one case that does need care, another bubble policy withdrawing a tie
    // this one still wants, is handled by handing the tie over (see
    // withdraw_ties()).
    if (up)
        return;

    build_ties(model);

}

template<typename TSeq>
inline void Bubbles<TSeq>::compute_partition(Model<TSeq> * model)
{

    bubble_id.assign(household_id.size(), -1);

    // Both rules read the contact network to decide which households may
    // bubble together. The ties a bubble policy has put there -- this one's,
    // or those of another policy on the same model -- are not contacts anybody
    // has: counting them would pull the next bubble towards the last one's
    // membership, and make one policy's grouping depend on another's. So
    // grouping skips them and sees the network the model actually has, which
    // is what keeps BubbleTies independent of BubbleFlavor.
    //
    // Skipping a tie rather than removing it keeps every other neighbor where
    // it was: add_edge() appends and rm_edge() preserves order, so the agents'
    // real contacts come up in the same order as if no tie had ever been
    // added, and the same draws group them the same way.
    HeldTies held(model, *this);

    if (flavor == BubbleFlavor::Household)
        partition_household(model, held);
    else
        partition_peer(model, held);

}

template<typename TSeq>
inline bool Bubbles<TSeq>::is_active(int today) const
{
    return (today >= start_day) && ((end_day < 0) || (today < end_day));
}

template<typename TSeq>
inline epiworld_double Bubbles<TSeq>::susceptibility_reduction(
    const Agent<TSeq> * p,
    const Agent<TSeq> * transmitter,
    Model<TSeq> * model
) const
{

    if (!is_active(static_cast< int >(model->today())))
        return 0.0;

    if (bubble_id.empty() || (p == nullptr) || (transmitter == nullptr))
        return 0.0;

    int bp = bubble_id[static_cast< size_t >(p->get_id())];
    int bt = bubble_id[static_cast< size_t >(transmitter->get_id())];

    if ((bp < 0) || (bt < 0))
        return 0.0;

    // Contacts inside the bubble are exactly what the policy keeps: they are
    // left alone.
    if (bp == bt)
        return 0.0;

    // Contacts outside the bubble are scaled by the transmission factor:
    // 0 = perfectly observed bubble (contact cut), 1 = the bubble imposes
    // nothing.
    epiworld_double factor = model->par(param_name);
    if (factor <= 0.0)
        return 1.0;
    if (factor >= 1.0)
        return 0.0;

    return static_cast<epiworld_double>(1.0) - factor;

}

template<typename TSeq>
inline void Bubbles<TSeq>::_setup(Model<TSeq> * model)
{

    if (household_id.size() != model->size())
        throw std::length_error(
            "Bubbles: household_id length (" +
            std::to_string(household_id.size()) +
            ") must equal the number of agents (" +
            std::to_string(model->size()) + ")."
        );

    if (ties == BubbleTies::Complete)
    {

        // Completing a bubble means editing ties at both ends, which is not a
        // meaningful thing to do to a directed network.
        if (model->is_directed())
            throw std::logic_error(
                "Bubbles: BubbleTies::Complete needs an undirected model."
            );

        // Rewiring moves ties *between* agents, so a tie this created could be
        // moved out from under it and would never be withdrawn. The two cannot
        // be combined. It is the function that matters, not the proportion:
        // Model::rewire() calls it on every step, and nothing obliges a
        // rewiring function to honour a proportion of zero.
        if (model->has_rewire_fun())
            throw std::logic_error(
                "Bubbles: BubbleTies::Complete cannot be combined with network "
                "rewiring; remove it with set_rewire_fun(nullptr)."
            );

    }

    // Whatever a previous run left in the network goes first. Normally there is
    // nothing to do -- a run withdraws its own ties on its last day (see
    // wants_ties) -- but a hand-driven day loop can stop early.
    restore_network(model);

    // ---- The transmission factor lives in the model -------------------------
    // The tool reads it on every exposure rather than holding a copy, so the
    // strictness of the policy can be inspected, calibrated, or switched
    // mid-run through the model's parameters. The value passed to the
    // constructor is only a default: a value already in the model (set by the
    // user, or read from a parameter file) is what governs the run.
    if (!model->has_param(param_name))
        model->add_param(transmission_factor, param_name);

    // ---- The partition ------------------------------------------------------
    // Drawn with the model's RNG, which the run has already seeded, so each
    // replicate of run_multiple() gets its own partition from its own seed and
    // never inherits one from a previous run.
    compute_partition(model);
    last_epoch = 0;

    // ---- The bubble as ties -------------------------------------------------
    // Under BubbleTies::Complete the bubble is not only a transmission rule: it
    // is a clique, so households are completed and merged households actually
    // meet. No queue bookkeeping is needed here -- Model::reset() has cleared
    // the queue and the dist_virus() events are still pending, so the
    // events_run() that follows this builds the counts from the finished
    // network.
    sync_ties(model);

    // ---- The tool: dampens out-of-bubble transmission -----------------------
    // It carries no state of its own: it finds the model's intervention by name
    // the first time it is used. It is registered without a distribution
    // function on purpose -- handing it out here, on every run, keeps the first
    // run and the ones after it (where the tool is already registered)
    // identical.
    if (!model->has_tool(this->get_name()))
    {
        BubbleTool<TSeq> bubble_tool(this->get_name(), this->get_name());
        model->add_tool(bubble_tool);
    }

    auto & bubble_tool = model->get_tool(this->get_name());
    for (size_t i = 0u; i < model->size(); ++i)
        model->get_agent(i).add_tool(*model, bubble_tool);

}

template<typename TSeq>
inline void Bubbles<TSeq>::reset(Model<TSeq> * model)
{

    // Model::reset() runs this once per run -- and once per replicate of
    // run_multiple(), on that replicate's own copy of the model -- just before
    // day 1, which is why the user has nothing to call: adding the intervention
    // to the model is the whole installation.
    this->model_id = static_cast< int >(model->get_sim_id());
    this->_setup(model);

}

template<typename TSeq>
inline void Bubbles<TSeq>::operator()(Model<TSeq> * model, int day)
{

    // Under Model::run() this never fires: reset() has already set us up for
    // this run. It is here for a model whose day loop is driven by hand, where
    // installing the policy a day late still beats running without it. The
    // simulation id is what tells one run from the next, so a copy of the model
    // (another replicate, another thread) sets itself up on its own.
    if (static_cast< int >(model->get_sim_id()) != this->model_id)
    {
        this->model_id = static_cast< int >(model->get_sim_id());
        this->_setup(model);
    }

    // Past setup, the daily event has two jobs: move the partition on at a
    // rewiring epoch, and keep the network matching the policy. The tool gates
    // itself by day. Because global events run after update_state(), both take
    // effect the following simulation step.
    if ((rewire_every > 0) && is_active(day))
    {

        int epoch = (day - start_day) / rewire_every;

        if (last_epoch != epoch)
        {

            // The standing clique belongs to the old partition, so it comes
            // down before the new one is drawn. Grouping would skip it anyway
            // (see compute_partition()), but walking ties that are about to go
            // only makes the grouping slower. Handing over rather than plainly
            // removing: another policy may want some of these ties for a
            // bubble of its own that is still open. sync_ties() below builds
            // the new partition's clique.
            withdraw_ties(model, true);

            compute_partition(model);
            last_epoch = epoch;

        }

    }

    // Cheap and immediate unless the policy is realized as ties: builds the
    // clique on the day the policy starts, redraws it when the partition moves
    // to a new epoch, and withdraws it when the policy lifts or the run ends.
    // Model::add_edge()/rm_edge() keep the queueing system in step, so this is
    // safe in the middle of a run.
    sync_ties(model);

}

template<typename TSeq>
inline Bubbles<TSeq> * Bubbles<TSeq>::get_from(
    Model<TSeq> & model,
    const std::string & name
)
{

    if (!model.has_globalevent(name))
        return nullptr;

    return dynamic_cast< Bubbles<TSeq> * >(&model.get_globalevent(name));

}

template<typename TSeq>
inline const std::vector< int > & Bubbles<TSeq>::get_bubble_id() const
{
    return bubble_id;
}

template<typename TSeq>
inline int Bubbles<TSeq>::get_last_epoch() const
{
    return last_epoch;
}

template<typename TSeq>
inline BubbleFlavor Bubbles<TSeq>::get_flavor() const
{
    return flavor;
}

template<typename TSeq>
inline BubbleTies Bubbles<TSeq>::get_ties() const
{
    return ties;
}

template<typename TSeq>
inline void Bubbles<TSeq>::set_ties(BubbleTies ties)
{
    this->ties = ties;
}

template<typename TSeq>
inline const std::vector< std::pair< size_t, size_t > > &
Bubbles<TSeq>::get_created_ties() const
{
    return created_ties;
}

template<typename TSeq>
inline const std::string & Bubbles<TSeq>::get_param_name() const
{
    return param_name;
}

template<typename TSeq>
inline std::unique_ptr< GlobalEvent<TSeq> > Bubbles<TSeq>::clone_ptr() const
{
    return std::make_unique< Bubbles<TSeq> >(*this);
}

#endif
