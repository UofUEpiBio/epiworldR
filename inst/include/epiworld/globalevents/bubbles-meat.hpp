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
    BubbleTies ties,
    std::string layer,
    double bubble_contact_prob
) :
    flavor(flavor),
    group_size(group_size),
    max_households(max_households),
    transmission_factor(transmission_factor),
    start_day(start_day),
    end_day(end_day),
    rewire_every(rewire_every),
    param_name(std::move(param_name)),
    ties(ties),
    layer(std::move(layer)),
    bubble_contact_prob(bubble_contact_prob)
{

    // Households by compact index, in order of first appearance
    std::unordered_map< size_t, size_t > index;
    household.reserve(household_id.size());
    for (size_t h : household_id)
        household.push_back(index.emplace(h, index.size()).first->second);

    n_households = index.size();

    this->set_name(name);
    this->set_day(-99); // runs at the end of every day

    if ((this->transmission_factor < 0.0) || (this->transmission_factor > 1.0))
        throw std::range_error(
            "Bubbles: transmission_factor must be in [0, 1]."
        );

    if (!((bubble_contact_prob >= 0.0) && (bubble_contact_prob <= 1.0)))
        throw std::range_error(
            "Bubbles: bubble_contact_prob must be in [0, 1]."
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

    size_t nh = n_households;

    // Build the household contact graph: households h1 and h2 are adjacent when
    // at least one member of h1 is connected to a member of h2 in the contact
    // network. Bubbles are grown along these ties.
    //
    // Grouping households that share NO tie would be a no-op: the intervention
    // can only suppress transmission along existing edges, never create new
    // ones, so bubbling two unconnected households changes nothing. Pairing at
    // random therefore degenerates to the household-only lockdown. It also
    // matches the policy being modelled: a household picks a bubble partner it
    // actually socialises with. Only network layers count: a bubble layer is
    // a mixing layer (see make_layer()).
    std::vector< std::vector< size_t > > hh_adj(nh);
    for (size_t a = 0u; a < household.size(); ++a)
    {
        size_t ha = household[a];
        for (auto * net : model->get_network_layers())
        for (auto * nb : net->neighbors_view(a))
        {
            size_t hb = household[static_cast< size_t >(nb->get_id())];
            if (hb != ha)
                hh_adj[ha].push_back(hb);
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
    for (size_t a = 0u; a < household.size(); ++a)
        bubble_id[a] = hh_bubble[household[a]];

}

template<typename TSeq>
inline void Bubbles<TSeq>::partition_peer(Model<TSeq> * model)
{

    size_t n  = household.size();
    size_t nh = n_households;

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

    std::vector< size_t > ext;

    for (size_t oi = 0u; oi < n; ++oi)
    {

        size_t a  = agent_order[oi];
        size_t ha = household[a];

        // This agent's household is already in a full bubble: it is out of the
        // pool and cannot take anyone else in.
        if (set_size[find(ha)] >= max_households)
            continue;

        // Households of this agent's contacts outside its own household
        ext.clear();
        for (auto * net : model->get_network_layers())
        for (auto * nb : net->neighbors_view(a))
        {
            size_t hb = household[static_cast< size_t >(nb->get_id())];
            if (hb != ha)
                ext.push_back(hb);
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
        size_t root = find(household[a]);
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
inline std::shared_ptr< layer::Mixing<TSeq> > Bubbles<TSeq>::make_layer(
    double household_contact_prob
) const
{

    if (layer.empty())
        throw std::logic_error(
            "Bubbles: make_layer() needs the layer's name; give it with "
            "set_ties(BubbleTies::Complete, name)."
        );

    // Everyone with their household: the bubbles are larger groups of the
    // same layer (households are never split)
    auto res = layer::mixing<TSeq>(layer);
    res->set_groups(household, n_households);
    res->set_contact_prob(household_contact_prob);

    return res;

}

template<typename TSeq>
inline layer::Mixing<TSeq> & Bubbles<TSeq>::bubble_layer(
    Model<TSeq> * model
) const
{

    layer::Mixing<TSeq> * res = (!layer.empty() && model->has_layer(layer)) ?
        model->get_layer(layer).as_mixing() : nullptr;

    if ((res == nullptr) || !res->within_groups())
        throw std::logic_error(
            "Bubbles: BubbleTies::Complete needs its within-group mixing layer "
            "'" + layer + "' in the model; add it with "
            "model.add_layer(bubbles.make_layer())."
        );

    return *res;

}

template<typename TSeq>
inline void Bubbles<TSeq>::sync_groups(Model<TSeq> * model)
{

    // Global events run after the day's transitions, so the moves take effect
    // at the end of the step (the layer applies them in its update): what
    // matters is whether the policy applies the next day, which is also when
    // the tool starts (or stops) damping.
    bool want = (ties == BubbleTies::Complete) &&
        is_active(static_cast< int >(model->today()) + 1);

    // The partition agents should be in (-1: their homes), and the one they are
    int target = want ? last_epoch : -1;
    if (target == placed_epoch)
        return;

    auto & groups = bubble_layer(model);
    for (size_t i = 0u; i < household.size(); ++i)
        groups.set_group(
            i, want ? static_cast< size_t >(bubble_id[i]) : household[i]
        );

    // One contact probability per group: a bubble that merged households
    // meets at the bubble's; a household on its own keeps its own
    if (want)
    {
        std::vector< double > q(home_prob.size(), bubble_contact_prob);
        std::vector< size_t > first(home_prob.size(), layer::Mixing<TSeq>::NO_GROUP);
        std::vector< bool > merged(home_prob.size(), false);
        for (size_t i = 0u; i < household.size(); ++i)
        {
            const size_t b = static_cast< size_t >(bubble_id[i]);
            if (first[b] == layer::Mixing<TSeq>::NO_GROUP)
                first[b] = household[i];
            else if (first[b] != household[i])
                merged[b] = true;
        }

        for (size_t b = 0u; b < q.size(); ++b)
            if (!merged[b] && (first[b] != layer::Mixing<TSeq>::NO_GROUP))
                q[b] = home_prob[first[b]];

        groups.set_contact_prob(std::move(q));
    }
    else
        groups.set_contact_prob(home_prob);

    placed_epoch = target;

}

template<typename TSeq>
inline void Bubbles<TSeq>::compute_partition(Model<TSeq> * model)
{

    bubble_id.assign(household.size(), -1);

    // Both rules read the network layers only, so a bubble layer (a mixing
    // layer) never pulls the next bubble towards the last one: BubbleTies
    // does not change the grouping.
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
inline void Bubbles<TSeq>::_setup(Model<TSeq> * model)
{

    if (household.size() != model->size())
        throw std::length_error(
            "Bubbles: household_id length (" +
            std::to_string(household.size()) +
            ") must equal the number of agents (" +
            std::to_string(model->size()) + ")."
        );

    if (ties == BubbleTies::Complete)
    {

        auto & groups = bubble_layer(model);

        // An agent is in one group per layer, so two policies need two layers
        for (size_t e = 0u; e < model->get_n_globalevents(); ++e)
        {
            const auto * other = dynamic_cast< const Bubbles<TSeq> * >(
                &model->get_globalevent(e)
            );

            if ((other != nullptr) && (other != this) &&
                (other->ties == BubbleTies::Complete) && (other->layer == layer))
                throw std::logic_error(
                    "Bubbles: two policies move agents in layer '" + layer +
                    "'; give each its own layer (make_layer())."
                );
        }

        // A move outside a run would edit the layer's baseline for good
        if (!model->is_running())
            throw std::logic_error(
                "Bubbles: BubbleTies::Complete runs within Model::run() (the "
                "layer's reset puts everyone back in their household)."
            );

        // The policy moves agents from their households and back, so the
        // layer must start from the household partition (as make_layer()
        // builds it)
        bool households = (groups.get_n_groups() == n_households);
        for (size_t i = 0u; households && (i < household.size()); ++i)
            households = (groups.get_group(i) == household[i]);

        if (!households)
            throw std::logic_error(
                "Bubbles: the layer '" + layer + "' does not have one group "
                "per household, as make_layer() builds it (" +
                std::to_string(groups.get_n_groups()) + " groups for " +
                std::to_string(n_households) + " households); add it with "
                "model.add_layer(bubbles.make_layer())."
            );

        // How often each household meets
        home_prob = groups.get_contact_prob();

    }

    placed_epoch = -1;

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

    // ---- The bubble in its layer -------------------------------------------
    // Under BubbleTies::Complete agents move to their bubble's group of the
    // layer; the moves take effect in the next() that follows the reset, so
    // the bubbles meet from day 1.
    sync_groups(model);

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
    // rewiring epoch, and keep the bubble layer matching the policy. The tool gates
    // itself by day. Because global events run after update_state(), both take
    // effect the following simulation step.
    if ((rewire_every > 0) && is_active(day))
    {

        int epoch = (day - start_day) / rewire_every;

        if (last_epoch != epoch)
        {

            compute_partition(model);
            last_epoch = epoch;

        }

    }

    // Moves agents into their bubbles on the day before the policy starts,
    // into the new ones at a rewiring epoch, and back when the policy lifts
    sync_groups(model);

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
inline void Bubbles<TSeq>::set_ties(
    BubbleTies ties,
    std::string layer,
    double bubble_contact_prob
)
{

    if (!((bubble_contact_prob >= 0.0) && (bubble_contact_prob <= 1.0)))
        throw std::range_error(
            "Bubbles: bubble_contact_prob must be in [0, 1]."
        );

    this->ties = ties;
    this->layer = std::move(layer);
    this->bubble_contact_prob = bubble_contact_prob;

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
