#ifndef EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP
#define EPIWORLD_GLOBALEVENTS_BUBBLES_MEAT_HPP

#include <algorithm>
#include <stdexcept>
#include <utility>
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

    // Shuffle household order (Fisher-Yates with the model RNG).
    std::vector< size_t > order(nh);
    for (size_t i = 0u; i < nh; ++i)
        order[i] = i;

    for (size_t i = nh; i > 1u; --i)
    {
        size_t j = static_cast< size_t >(model->runif_index(static_cast<uint32_t>(i)));
        std::swap(order[i - 1u], order[j]);
    }

    // Chunk consecutive households into bubbles of `group_size`.
    std::vector< int > hh_bubble(nh, -1);
    for (size_t pos = 0u; pos < nh; ++pos)
        hh_bubble[order[pos]] = static_cast< int >(pos / group_size);

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

    // Eagerly compute the initial partition so the tool is effective from the
    // first simulation step (the scheduler runs after agent updates each day).
    compute_partition(&model);
    state->last_sim_id = -1;
    state->last_epoch  = -1;

    // ---- The bubble tool: blocks cross-bubble transmission -----------------
    auto st  = state;
    int  sd  = start_day;
    int  ed  = end_day;
    epiworld_double wf = within_factor;

    Tool<TSeq> bubble_tool(name, 1.0, true);
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

    model.add_tool(bubble_tool);

    // ---- The scheduler: (re)computes the partition per replicate/epoch -----
    Bubbles<TSeq> self = *this; // shares `state` via the shared_ptr
    model.add_globalevent(
        [self](Model<TSeq> * m) -> void {

            int today = static_cast< int >(m->today());

            bool on = (today >= self.start_day) &&
                ((self.end_day < 0) || (today < self.end_day));
            if (!on)
                return;

            int sim = static_cast< int >(m->get_sim_id());
            int epoch = (self.rewire_every > 0) ?
                ((today - self.start_day) / self.rewire_every) : 0;

            // Already up to date for this (replicate, epoch).
            if ((self.state->last_sim_id == sim) &&
                (self.state->last_epoch == epoch))
                return;

            // First run after deploy already holds the epoch-0 partition; claim
            // it without recomputing so a fixed policy stays reproducible.
            if ((self.state->last_sim_id == -1) && (epoch == 0))
            {
                self.state->last_sim_id = sim;
                self.state->last_epoch  = epoch;
                return;
            }

            self.compute_partition(m);
            self.state->last_sim_id = sim;
            self.state->last_epoch  = epoch;

        },
        name + " (scheduler)",
        -99
    );

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
