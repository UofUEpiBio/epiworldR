#ifndef EPIWORLD_PERSON_MEAT_HPP
#define EPIWORLD_PERSON_MEAT_HPP

// template<typename Ta>
// inline bool IN(Ta & a, std::vector< Ta > & b);

#define CHECK_COALESCE_(proposed_, virus_tool_, alt_) \
    if (static_cast<int>(proposed_) == -99) {\
        if (static_cast<int>(virus_tool_) == -99) \
            (proposed_) = (alt_);\
        else (proposed_) = (virus_tool_);}

// To large to add directly here
#include "agent-actions-meat.hpp"

template<typename TSeq>
inline Agent<TSeq>::Agent() {}

template<typename TSeq>
inline Agent<TSeq>::Agent(Agent<TSeq> && p) :
    model(p.model),
    neighbors(std::move(p.neighbors)),
    neighbors_locations(std::move(p.neighbors_locations)),
    n_neighbors(p.n_neighbors),
    entities(std::move(p.entities)),
    entities_locations(std::move(p.entities_locations)),
    n_entities(p.n_entities),
    state(p.state),
    state_prev(p.state_prev), 
    state_last_changed(p.state_last_changed),
    id(p.id),
    viruses(std::move(p.viruses)),  /// Needs to be adjusted
    n_viruses(p.n_viruses),
    tools(std::move(p.tools)), /// Needs to be adjusted
    n_tools(p.n_tools),
    add_virus_(std::move(p.add_virus_)),
    add_tool_(std::move(p.add_tool_)),
    add_entity_(std::move(p.add_entity_)),
    rm_virus_(std::move(p.rm_virus_)),
    rm_tool_(std::move(p.rm_tool_)),
    rm_entity_(std::move(p.rm_entity_)),
    action_counter(p.action_counter)
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus

    int loc = 0;
    for (auto & v : viruses)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        v->agent = this;
        v->pos_in_agent = loc++;

    }

    loc = 0;
    for (auto & t : tools)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        t->agent     = this;
        t->pos_in_agent = loc++;

    }
    
}

template<typename TSeq>
inline Agent<TSeq>::Agent(const Agent<TSeq> & p) :
    model(p.model),
    neighbors(p.neighbors),
    neighbors_locations(p.neighbors_locations),
    n_neighbors(p.n_neighbors),
    entities(p.entities),
    entities_locations(p.entities_locations),
    n_entities(p.n_entities),
    sampled_agents(0u),
    sampled_agents_n(0u),
    sampled_agents_left_n(0u),
    date_last_build_sample(-99)
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus
    viruses.reserve(p.get_n_viruses());
    n_viruses = viruses.size();
    for (size_t i = 0u; i < n_viruses; ++i)
    {

        // Will create a copy of the virus, with the exeption of
        // the virus code
        viruses.emplace_back(std::make_shared<Virus<TSeq>>(*p.viruses[i]));
        viruses.back()->set_agent(this, i);

    }

    tools.reserve(p.get_n_tools());
    n_tools = tools.size();
    for (size_t i = 0u; i < n_tools; ++i)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        tools.emplace_back(std::make_shared<Tool<TSeq>>(*p.tools[i]));
        tools.back()->set_agent(this, i);

    }

    add_virus_ = p.add_virus_;
    add_tool_  = p.add_tool_;
    rm_virus_  = p.rm_virus_;
    rm_tool_   = p.rm_tool_;
    
}

template<typename TSeq>
inline Agent<TSeq> & Agent<TSeq>::operator=(
    const Agent<TSeq> & other_agent
) 
{

    model = other_agent.model;

    neighbors = other_agent.neighbors;
    neighbors_locations = other_agent.neighbors_locations;
    n_neighbors = other_agent.n_neighbors;

    entities = other_agent.entities;
    entities_locations = other_agent.entities_locations;
    n_entities = other_agent.n_entities;

    sampled_agents.clear();
    sampled_agents_n = 0;
    sampled_agents_left_n = 0;
    date_last_build_sample = -99;

    // neighbors           = other_agent.neighbors;
    // entities            = other_agent.entities;
    // entities_locations  = other_agent.entities_locations;
    // n_entities          = other_agent.n_entities;
    state              = other_agent.state;
    state_prev         = other_agent.state_prev;
    state_last_changed = other_agent.state_last_changed;
    id                  = other_agent.id;
    
    // viruses             = other_agent.viruses;
    n_viruses           = other_agent.n_viruses;
    viruses.resize(n_viruses, nullptr);
    for (size_t i = 0u; i < n_viruses; ++i)
    {
        viruses[i] = std::make_shared<Virus<TSeq>>(*other_agent.viruses[i]);
        viruses[i]->set_agent(this, i);
    }
    
    // tools               = other_agent.tools;
    n_tools             = other_agent.n_tools;
    for (size_t i = 0u; i < n_tools; ++i)
    {
        tools[i] = std::make_shared<Tool<TSeq>>(*other_agent.tools[i]);
        tools[i]->set_agent(this, i);
    }

    add_virus_          = other_agent.add_virus_;
    add_tool_           = other_agent.add_tool_;
    add_entity_         = other_agent.add_entity_;
    rm_virus_           = other_agent.rm_virus_;
    rm_tool_            = other_agent.rm_tool_;
    rm_entity_          = other_agent.rm_entity_;
    action_counter      = other_agent.action_counter;
    
    return *this;
    
}

template<typename TSeq>
inline void Agent<TSeq>::add_tool(
    ToolPtr<TSeq> tool,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
) {

    // Checking the virus exists
    if (tool->get_id() >= static_cast<int>(model->get_db().get_n_tools()))
        throw std::range_error("The tool with id: " + std::to_string(tool->get_id()) + 
            " has not been registered. There are only " + std::to_string(model->get_n_tools()) + 
            " included in the model.");
    

    model->actions_add(
        this, nullptr, tool, nullptr, state_new, queue, add_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::add_tool(
    Tool<TSeq> tool,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    ToolPtr<TSeq> tool_ptr = std::make_shared< Tool<TSeq> >(tool);
    add_tool(tool_ptr, model, state_new, queue);
}

template<typename TSeq>
inline void Agent<TSeq>::add_virus(
    VirusPtr<TSeq> virus,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    // Checking the virus exists
    if (virus->get_id() >= static_cast<int>(model->get_db().get_n_viruses()))
        throw std::range_error("The virus with id: " + std::to_string(virus->get_id()) + 
            " has not been registered. There are only " + std::to_string(model->get_n_viruses()) + 
            " included in the model.");

    model->actions_add(
        this, virus, nullptr, nullptr, state_new, queue, add_virus_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::add_virus(
    Virus<TSeq> virus,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    VirusPtr<TSeq> virus_ptr = std::make_shared< Virus<TSeq> >(virus);
    add_virus(virus_ptr, model, state_new, queue);
}

template<typename TSeq>
inline void Agent<TSeq>::add_entity(
    Entity<TSeq> & entity,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (model != nullptr)
    {

        model->actions_add(
            this, nullptr, nullptr, &entity, state_new, queue, add_entity_, -1, -1
        );

    }
    else // If no model is passed, then we assume that we only need to add the
         // model entity
    {

        Action<TSeq> a(
                this, nullptr, nullptr, &entity, state_new, queue, add_entity_,
                -1, -1
            );

        default_add_entity(a, model); /* passing model makes nothing */

    }

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    epiworld_fast_uint tool_idx,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool_idx >= n_tools)
        throw std::range_error(
            "The Tool you want to remove is out of range. This Agent only has " +
            std::to_string(n_tools) + " tools."
        );

    model->actions_add(
        this, nullptr, tools[tool_idx], nullptr, state_new, queue, rm_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    ToolPtr<TSeq> & tool,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model->actions_add(
        this, nullptr, tool, nullptr, state_new, queue, rm_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    epiworld_fast_uint virus_idx,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    if (virus_idx >= n_viruses)
        throw std::range_error(
            "The Virus you want to remove is out of range. This Agent only has " +
            std::to_string(n_viruses) + " viruses."
        );
    else if (n_viruses == 0u)
        throw std::logic_error(
            "There is no virus to remove here!"
        );

    #ifdef EPI_DEBUG
    if (viruses[virus_idx]->pos_in_agent >= static_cast<int>(n_viruses))
    {
        throw std::logic_error(
            "[epi-debug]::rm_virus the position in the agent is wrong."
            );
    }
    #endif

    model->actions_add(
        this, viruses[virus_idx], nullptr, nullptr, state_new, queue,
        default_rm_virus<TSeq>, -1, -1
        );
    
}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    VirusPtr<TSeq> & virus,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (virus->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model->actions_add(
        this, virus, nullptr, nullptr, state_new, queue,
        default_rm_virus<TSeq>, -1, -1
        );


}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    epiworld_fast_uint entity_idx,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (entity_idx >= n_entities)
        throw std::range_error(
            "The Entity you want to remove is out of range. This Agent only has " +
            std::to_string(n_entities) + " entitites."
        );
    else if (n_entities == 0u)
        throw std::logic_error(
            "There is entity to remove here!"
        );

    model->actions_add(
        this, nullptr, nullptr, model->entities[entity_idx], state_new, queue, 
        default_rm_entity, entities_locations[entity_idx], entity_idx
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    Entity<TSeq> & entity,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    // Looking for entity location in the agent
    int entity_idx = -1;
    for (size_t i = 0u; i < n_entities; ++i)
    {
        if (entities[i] == entity->get_id())
            entity_idx = i;
    }

    if (entity_idx == -1)
        throw std::logic_error(
            "The agent " + std::to_string(id) + " is not associated with entity \"" +
            entity.get_name() + "\"."
            );


    model->actions_add(
        this, nullptr, nullptr, entities[entity_idx], state_new, queue, 
        default_rm_entity, entities_locations[entity_idx], entity_idx
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_agent_by_virus(
    epiworld_fast_uint virus_idx,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (state_new == -99)
        state_new = state;

    if (virus_idx >= n_viruses)
        throw std::range_error(
            std::string("The virus trying to remove the agent is out of range. ") +
            std::string("This agent has only ") + std::to_string(n_viruses) + 
            std::string(" and you are trying to remove virus # ") +
            std::to_string(virus_idx) + std::string(".")
            );

    // Removing viruses
    for (size_t i = 0u; i < n_viruses; ++i)
    {
        if (i != virus_idx)
            rm_virus(i, model);
    }

    // Changing state to new_state
    VirusPtr<TSeq> & v = viruses[virus_idx];
    epiworld_fast_int dead_state, dead_queue;
    v->get_state(nullptr, nullptr, &dead_state);
    v->get_queue(nullptr, nullptr, &dead_queue);

    if (queue != -99)
        dead_queue = queue;

    change_state(
        model,
        // Either preserve the current state or apply a new one
        (dead_state < 0) ? state : static_cast<epiworld_fast_uint>(dead_state),

        // By default, it will be removed from the queue... unless the user
        // says the contrary!
        (dead_queue == -99) ? Queue<TSeq>::NoOne : dead_queue
    );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_agent_by_virus(
    VirusPtr<TSeq> & virus,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (virus->get_agent() == nullptr)
        throw std::logic_error("The virus trying to remove the agent has no host.");

    if (virus->get_agent()->id != id)
        throw std::logic_error("Viruses can only remove their hosts'.");

    rm_agent_by_virus(
        virus->pos_in_agent,
        model,
        state_new,
        queue
    );

}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {

    return model->susceptibility_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->transmission_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->recovery_enhancer_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_death_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->death_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline int Agent<TSeq>::get_id() const
{
    return id;
}

template<typename TSeq>
inline Viruses<TSeq> Agent<TSeq>::get_viruses() {

    return Viruses<TSeq>(*this);

}

template<typename TSeq>
inline const Viruses_const<TSeq> Agent<TSeq>::get_viruses() const {

    return Viruses_const<TSeq>(*this);
    
}

template<typename TSeq>
inline VirusPtr<TSeq> & Agent<TSeq>::get_virus(int i) {
    return viruses.at(i);
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_viruses() const noexcept
{
    return n_viruses;
}

template<typename TSeq>
inline Tools<TSeq> Agent<TSeq>::get_tools() {
    return Tools<TSeq>(*this);
}

template<typename TSeq>
inline const Tools_const<TSeq> Agent<TSeq>::get_tools() const {
    return Tools_const<TSeq>(*this);
}

template<typename TSeq>
inline ToolPtr<TSeq> & Agent<TSeq>::get_tool(int i)
{
    return tools.at(i);
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_tools() const noexcept
{
    return n_tools;
}

template<typename TSeq>
inline void Agent<TSeq>::mutate_virus()
{

    for (auto & v : viruses)
        v->mutate();

}

template<typename TSeq>
inline void Agent<TSeq>::add_neighbor(
    Agent<TSeq> & p,
    bool check_source,
    bool check_target
) {
    // Can we find the neighbor?
    bool found = false;
    if (check_source)
    {

        for (auto & n: neighbors)    
            if (static_cast<int>(n) == p.get_id())
            {
                found = true;
                break;
            }

    }

    // Three things going on here:
    // - Where in the neighbor will this be
    // - What is the neighbor's id
    // - Increasing the number of neighbors
    if (!found)
    {

        neighbors_locations.push_back(p.get_n_neighbors());
        neighbors.push_back(p.get_id());
        n_neighbors++;

    }


    found = false;
    if (check_target)
    {

        for (auto & n: p.neighbors)
            if (static_cast<int>(n) == id)
            {
                found = true;
                break;
            }
    
    }

    if (!found)
    {

        p.neighbors_locations.push_back(n_neighbors - 1);
        p.neighbors.push_back(id);
        p.n_neighbors++;
        
    }
    

}

template<typename TSeq>
inline void Agent<TSeq>::swap_neighbors(
    Agent<TSeq> & other,
    size_t n_this,
    size_t n_other
)
{

    // Getting the agents
    auto & pop = model->population;
    auto & neigh_this  = pop[neighbors[n_this]];
    auto & neigh_other = pop[other.neighbors[n_other]];

    // Getting the locations in the neighbors
    size_t loc_this_in_neigh = neighbors_locations[n_this];
    size_t loc_other_in_neigh = other.neighbors_locations[n_other];

    // Changing ids
    std::swap(neighbors[n_this], other.neighbors[n_other]);

    if (!model->directed)
    {
        std::swap(
            neigh_this.neighbors[loc_this_in_neigh],
            neigh_other.neighbors[loc_other_in_neigh]
            );

        // Changing the locations
        std::swap(neighbors_locations[n_this], other.neighbors_locations[n_other]);
        
        std::swap(
            neigh_this.neighbors_locations[loc_this_in_neigh],
            neigh_other.neighbors_locations[loc_other_in_neigh]
            );
    }

}

template<typename TSeq>
inline std::vector< Agent<TSeq> *> Agent<TSeq>::get_neighbors()
{
    std::vector< Agent<TSeq> * > res(n_neighbors, nullptr);
    for (size_t i = 0u; i < n_neighbors; ++i)
        res[i] = &model->population[neighbors[i]];

    return res;
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_neighbors() const
{
    return n_neighbors;
}

template<typename TSeq>
inline void Agent<TSeq>::change_state(
    Model<TSeq> * model,
    epiworld_fast_uint new_state,
    epiworld_fast_int queue
    )
{

    model->actions_add(
        this, nullptr, nullptr, nullptr, new_state, queue, nullptr, -1, -1
    );
    
    return;

}

template<typename TSeq>
inline const epiworld_fast_uint & Agent<TSeq>::get_state() const {
    return state;
}

template<typename TSeq>
inline void Agent<TSeq>::reset()
{

    this->viruses.clear();
    n_viruses = 0u;

    this->tools.clear();
    n_tools = 0u;

    this->state = 0u;
    this->state_prev = 0u;

    this->state_last_changed = -1;
    
}

template<typename TSeq>
inline bool Agent<TSeq>::has_tool(epiworld_fast_uint t) const
{

    for (auto & tool : tools)
        if (tool->get_id() == static_cast<int>(t))
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_tool(std::string name) const
{

    for (auto & tool : tools)
        if (tool->get_name() == name)
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_tool(const Tool<TSeq> & tool) const
{

    return has_tool(tool.get_id());

}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(epiworld_fast_uint t) const
{
    for (auto & v : viruses)
        if (v->get_id() == static_cast<int>(t))
            return true;

    return false;
}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(std::string name) const
{
    
    for (auto & v : viruses)
        if (v->get_name() == name)
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(const Virus<TSeq> & virus) const
{

    return has_virus(virus.get_id());

}

template<typename TSeq>
inline void Agent<TSeq>::print(
    Model<TSeq> * model,
    bool compressed
    ) const
{

    if (compressed)
    {
        printf_epiworld(
            "Agent: %i, state: %s (%lu), Nvirus: %lu, NTools: %lu, NNeigh: %lu\n",
            id, model->states_labels[state].c_str(), state, n_viruses, n_tools, neighbors.size()
        );
    }
    else {

        printf_epiworld("Information about agent id %i\n", this->id);
        printf_epiworld("  State        : %s (%lu)\n", model->states_labels[state].c_str(), state);
        printf_epiworld("  Virus count  : %lu\n", n_viruses);
        printf_epiworld("  Tool count   : %lu\n", n_tools);
        printf_epiworld("  Neigh. count : %lu\n", neighbors.size());

        size_t nfeats = model->get_agents_data_ncols();
        if (nfeats > 0)
        {

            printf_epiworld("This model includes features (%lu): [ ", nfeats);

            int max_to_show = static_cast<int>((nfeats > 10)? 10 : nfeats);

            for (int k = 0; k < max_to_show; ++k)
            {
                printf_epiworld("%.2f", this->operator[](k));

                if (k != (max_to_show - 1))
                {
                    printf_epiworld(", ");
                } else {
                    printf_epiworld(" ]\n");
                }

            }
            
        }

    }

    return;

}

template<typename TSeq>
inline double & Agent<TSeq>::operator()(size_t j)
{

    if (model->agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model->agents_data + j * model->size() + id);

}

template<typename TSeq>
inline double & Agent<TSeq>::operator[](size_t j)
{
    return *(model->agents_data + j * model->size() + id);
}

template<typename TSeq>
inline double Agent<TSeq>::operator()(size_t j) const
{

    if (model->agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model->agents_data + j * model->size() + id);

}

template<typename TSeq>
inline double Agent<TSeq>::operator[](size_t j) const
{
    return *(model->agents_data + j * model->size() + id);
}

template<typename TSeq>
inline Entities<TSeq> Agent<TSeq>::get_entities()
{
    return Entities<TSeq>(*this);
}

template<typename TSeq>
inline const Entities_const<TSeq> Agent<TSeq>::get_entities() const
{
    return Entities_const<TSeq>(*this);
}

template<typename TSeq>
inline const Entity<TSeq> & Agent<TSeq>::get_entity(size_t i) const
{
    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->entities[entities[i]];
}

template<typename TSeq>
inline Entity<TSeq> & Agent<TSeq>::get_entity(size_t i)
{
    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->entities[entities[i]];
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_entities() const
{
    return n_entities;
}

template<typename TSeq>
inline bool Agent<TSeq>::operator==(const Agent<TSeq> & other) const
{

    EPI_DEBUG_FAIL_AT_TRUE(
        n_neighbors != other.n_neighbors,
        "Agent:: n_eighbors don't match"
        )

    
    for (size_t i = 0u; i < n_neighbors; ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            neighbors[i] != other.neighbors[i],
            "Agent:: neighbor[i] don't match"
        )
    }
    
    EPI_DEBUG_FAIL_AT_TRUE(
        n_entities != other.n_entities,
        "Agent:: n_entities don't match"
        )
    
    
    for (size_t i = 0u; i < n_entities; ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            entities[i] != other.entities[i],
            "Agent:: entities[i] don't match"
        )
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        state != other.state,
        "Agent:: state don't match"
        )
        

    EPI_DEBUG_FAIL_AT_TRUE(
        state_prev != other.state_prev,
        "Agent:: state_prev don't match"
        )
        

    // EPI_DEBUG_FAIL_AT_TRUE(
    //     state_last_changed != other.state_last_changed,
    //     "Agent:: state_last_changed don't match"
    //     ) ///< Last time the agent was updated.
        
    
    EPI_DEBUG_FAIL_AT_TRUE(
        n_viruses != other.n_viruses,
        "Agent:: n_viruses don't match"
        )
        

    for (size_t i = 0u; i < n_viruses; ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            *viruses[i] != *other.viruses[i],
            "Agent:: viruses[i] don't match"
        )
         
    }

    EPI_DEBUG_FAIL_AT_TRUE(n_tools != other.n_tools, "Agent:: n_tools don't match")

    for (size_t i = 0u; i < n_tools; ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            tools[i] != other.tools[i],
            "Agent:: tools[i] don't match"
        )
         
    }   
    
    return true;
    
}

#undef CHECK_COALESCE_

#endif