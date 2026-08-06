#ifndef EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP
#define EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP

// Standard library headers are included at global scope by epiworld.hpp (this
// file is included from within `namespace epiworld`, so system headers must not
// be re-included here).
#include "bubbles-bones.hpp"

template<typename TSeq>
inline Bubbles<TSeq>::Bubbles(
    std::vector< size_t > household_id,
    BubbleFlavor flavor,
    size_t group_size,
    epiworld_double within_factor,
    int start_day,
    int end_day,
    int rewire_every,
    std::string name
) :
    household_id(std::move(household_id)),
    flavor(flavor),
    group_size(group_size),
    within_factor(within_factor),
    start_day(start_day),
    end_day(end_day),
    rewire_every(rewire_every),
    name(std::move(name)),
    state(std::make_shared< BubbleState >())
{

    if ((this->within_factor < 0.0) || (this->within_factor > 1.0))
        throw std::range_error(
            "Bubbles: within_factor must be in [0, 1]."
        );

    if ((flavor == BubbleFlavor::Household) && (group_size < 1u))
        throw std::range_error(
            "Bubbles: group_size (households per bubble) must be >= 1."
        );

    if ((this->end_day >= 0) && (this->end_day <= this->start_day))
        throw std::range_error(
            "Bubbles: end_day must be greater than start_day (or negative)."
        );

}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_household(Model<TSeq> * model) const
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
    std::vector< std::vector< size_t > > hh_adj(nh);
    auto & pop_all = model->get_agents();
    for (size_t a = 0u; a < household_id.size(); ++a)
    {
        size_t ha = hh_index[household_id[a]];
        for (auto * nb : pop_all[a].get_neighbors(*model))
        {
            size_t b = static_cast< size_t >(nb->get_id());
            if (household_id[a] == household_id[b])
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
        state->bubble_id[a] = hh_bubble[hh_index[household_id[a]]];

}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_peer(Model<TSeq> * model) const
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

    // Disjoint-set (union-find) over households.
    std::vector< size_t > parent(nh);
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

    // Each agent selects up to `group_size` distinct external peers (neighbors
    // in a different household); each selection merges the two households.
    auto & pop = model->get_agents();
    for (size_t a = 0u; a < n; ++a)
    {

        std::vector< size_t > ext;
        for (auto * nb : pop[a].get_neighbors(*model))
        {
            size_t nid = static_cast< size_t >(nb->get_id());
            if (household_id[nid] != household_id[a])
                ext.push_back(nid);
        }

        size_t k = std::min(group_size, ext.size());
        for (size_t i = 0u; i < k; ++i)
        {
            // Partial Fisher-Yates: draw the i-th distinct peer.
            size_t j = i + static_cast< size_t >(
                model->runif_index(static_cast<uint32_t>(ext.size() - i))
            );
            std::swap(ext[i], ext[j]);

            size_t ra = find(hh_index[household_id[a]]);
            size_t rb = find(hh_index[household_id[ext[i]]]);
            if (ra != rb)
                parent[ra] = rb;
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
            state->bubble_id[a] = next_label;
            ++next_label;
        }
        else
        {
            state->bubble_id[a] = it->second;
        }
    }

}

template<typename TSeq>
inline void Bubbles<TSeq>::compute_partition(Model<TSeq> * model) const
{

    state->bubble_id.assign(household_id.size(), -1);

    if (flavor == BubbleFlavor::Household)
        partition_household(model);
    else
        partition_peer(model);

}

template<typename TSeq>
inline void Bubbles<TSeq>::deploy(Model<TSeq> & model)
{

    if (household_id.size() != model.size())
        throw std::length_error(
            "Bubbles: household_id length (" +
            std::to_string(household_id.size()) +
            ") must equal the number of agents (" +
            std::to_string(model.size()) + ")."
        );

    // ---- The bubble tool: blocks cross-bubble transmission -----------------
    auto st  = state;
    int  sd  = start_day;
    int  ed  = end_day;
    epiworld_double wf = within_factor;

    Tool<TSeq> bubble_tool(name);
    bubble_tool.set_susceptibility_reduction_fun(
        [st, sd, ed, wf](
            Tool<TSeq> &,
            Agent<TSeq> * p,
            VirusPtr<TSeq> & v,
            Model<TSeq> * m
        ) -> epiworld_double {

            int today = static_cast< int >(m->today());

            // Policy window.
            if (today < sd)
                return 0.0;
            if ((ed >= 0) && (today >= ed))
                return 0.0;

            if (st->bubble_id.empty())
                return 0.0;

            Agent<TSeq> * transmitter = v->get_agent();
            if (transmitter == nullptr)
                return 0.0;

            int bp = st->bubble_id[static_cast< size_t >(p->get_id())];
            int bt = st->bubble_id[static_cast< size_t >(transmitter->get_id())];

            if ((bp < 0) || (bt < 0))
                return 0.0;

            // Same bubble: within-bubble distancing. Different bubble: block.
            return (bp == bt) ?
                (static_cast<epiworld_double>(1.0) - wf) :
                static_cast<epiworld_double>(1.0);

        }
    );

    // ---- Distribution: recompute the partition at reset time ---------------
    // The tool's distribution function runs from Model::reset() -> dist_tools(),
    // i.e. *before* day 1 and with the run's (per-replicate) RNG already seeded.
    // Computing the partition here — rather than eagerly in deploy() — keeps the
    // epoch-0 partition fresh for every replicate of run_multiple() and avoids
    // any dependence on a stale partition left over from a previous run. It then
    // distributes the tool to every agent (prevalence 1.0).
    Bubbles<TSeq> self = *this; // shares `state` via the shared_ptr
    ToolToAgentFun<TSeq> distribute_all = distribute_tool_randomly<TSeq>(1.0, true);
    bubble_tool.set_distribution(
        [self, distribute_all](Tool<TSeq> & tool, Model<TSeq> * m) -> void {
            self.compute_partition(m);
            self.state->last_sim_id = static_cast< int >(m->get_sim_id());
            self.state->last_epoch  = 0;
            distribute_all(tool, m);
        }
    );

    model.add_tool(bubble_tool);

    // ---- The scheduler: re-randomizes the partition at rewiring epochs -----
    // The epoch-0 partition is installed at reset (above); this event only
    // handles rewire_every > 0, recomputing when the epoch advances. Because
    // global events run after update_state(), a rewired partition takes effect
    // the following simulation step. Deactivation (end_day) needs no event: the
    // tool gates itself by day.
    if (rewire_every > 0)
    {
        model.add_globalevent(
            [self](Model<TSeq> * m) -> void {

                int today = static_cast< int >(m->today());

                bool on = (today >= self.start_day) &&
                    ((self.end_day < 0) || (today < self.end_day));
                if (!on)
                    return;

                int sim   = static_cast< int >(m->get_sim_id());
                int epoch = (today - self.start_day) / self.rewire_every;

                // Already up to date for this (replicate, epoch).
                if ((self.state->last_sim_id == sim) &&
                    (self.state->last_epoch == epoch))
                    return;

                self.compute_partition(m);
                self.state->last_sim_id = sim;
                self.state->last_epoch  = epoch;

            },
            name + " (scheduler)",
            -99
        );
    }

}

template<typename TSeq>
inline std::shared_ptr< BubbleState > Bubbles<TSeq>::get_state() const
{
    return state;
}

template<typename TSeq>
inline const std::vector< int > & Bubbles<TSeq>::get_bubble_id() const
{
    return state->bubble_id;
}

template<typename TSeq>
inline BubbleFlavor Bubbles<TSeq>::get_flavor() const
{
    return flavor;
}

#endif
