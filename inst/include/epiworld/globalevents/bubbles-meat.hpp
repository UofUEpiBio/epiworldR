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
    std::string param_name
) :
    household_id(std::move(household_id)),
    flavor(flavor),
    group_size(group_size),
    max_households(max_households),
    transmission_factor(transmission_factor),
    start_day(start_day),
    end_day(end_day),
    rewire_every(rewire_every),
    param_name(std::move(param_name))
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
inline void Bubbles<TSeq>::partition_household(Model<TSeq> * model)
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
        bubble_id[a] = hh_bubble[hh_index[household_id[a]]];

}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_peer(Model<TSeq> * model)
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

        // Households of this agent's contacts outside its own household.
        ext.clear();
        for (auto * nb : pop[a].get_neighbors(*model))
        {
            size_t nid = static_cast< size_t >(nb->get_id());
            if (household_id[nid] != household_id[a])
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
inline void Bubbles<TSeq>::compute_partition(Model<TSeq> * model)
{

    bubble_id.assign(household_id.size(), -1);

    if (flavor == BubbleFlavor::Household)
        partition_household(model);
    else
        partition_peer(model);

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
inline void Bubbles<TSeq>::operator()(Model<TSeq> * model, int day)
{

    // Only rewiring needs a daily event; the tool gates itself by day, and the
    // epoch-0 partition is computed at reset (see deploy()). Because global
    // events run after update_state(), a rewired partition takes effect the
    // following simulation step.
    if (rewire_every <= 0)
        return;

    if (!is_active(day))
        return;

    int sim   = static_cast< int >(model->get_sim_id());
    int epoch = (day - start_day) / rewire_every;

    // Already up to date for this (replicate, epoch).
    if ((last_sim_id == sim) && (last_epoch == epoch))
        return;

    compute_partition(model);
    last_sim_id = sim;
    last_epoch  = epoch;

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

    // ---- The transmission factor lives in the model -------------------------
    // The tool reads it on every exposure rather than holding a copy, so the
    // strictness of the policy can be inspected, calibrated, or switched
    // mid-run through the model's parameters.
    model.add_param(transmission_factor, param_name, true);

    // ---- The intervention itself -------------------------------------------
    // add_globalevent() clones this object, so the model owns the partition and
    // the rewiring schedule. Copies of the model (e.g. the per-thread copies
    // made by run_multiple) clone it again and are fully independent.
    model.add_globalevent(*this);

    // ---- The tool: dampens out-of-bubble transmission -----------------------
    // It carries no state of its own: it finds the model's intervention by name
    // the first time it is used.
    BubbleTool<TSeq> bubble_tool(this->get_name(), this->get_name());

    // ---- Distribution: recompute the partition at reset time ----------------
    // The tool's distribution function runs from Model::reset() -> dist_tools(),
    // i.e. *before* day 1 and with the run's (per-replicate) RNG already seeded.
    // Computing the partition here — rather than eagerly in deploy() — keeps the
    // epoch-0 partition fresh for every replicate of run_multiple() and avoids
    // any dependence on a stale partition left over from a previous run. It then
    // distributes the tool to every agent (prevalence 1.0).
    std::string ename = this->get_name();
    ToolToAgentFun<TSeq> distribute_all = distribute_tool_randomly<TSeq>(1.0, true);
    bubble_tool.set_distribution(
        [ename, distribute_all](Tool<TSeq> & tool, Model<TSeq> * m) -> void {

            Bubbles<TSeq> * self = Bubbles<TSeq>::get_from(*m, ename);

            if (self == nullptr)
                throw std::logic_error(
                    "Bubbles: the intervention '" + ename +
                    "' is not installed on this model."
                );

            self->compute_partition(m);
            self->last_sim_id = static_cast< int >(m->get_sim_id());
            self->last_epoch  = 0;

            distribute_all(tool, m);

        }
    );

    model.add_tool(bubble_tool);

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
