#ifndef EPIWORLD_PERSON_MEAT_HPP
#define EPIWORLD_PERSON_MEAT_HPP

#define CHECK_COALESCE_(proposed_, virus_tool_, alt_) \
    if (static_cast<int>(proposed_) == -99) {\
        if (static_cast<int>(virus_tool_) == -99) \
            (proposed_) = (alt_);\
        else (proposed_) = (virus_tool_);}

// To large to add directly here
#include "agent-events-meat.hpp"

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
    tools(std::move(p.tools)), /// Needs to be adjusted
    n_tools(p.n_tools)
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus
    if (p.virus != nullptr)
    {
        virus = std::move(p.virus);
        virus->set_agent(this);
    }

    int loc = 0;
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
    n_neighbors(p.n_neighbors),
    entities(p.entities),
    entities_locations(p.entities_locations),
    n_entities(p.n_entities)
{

    if (n_neighbors > 0u)
    {
        neighbors = new std::vector< size_t >(*p.neighbors);
        neighbors_locations = new std::vector< size_t >(*p.neighbors_locations);
    }

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus
    if (p.virus != nullptr)
    {
        virus = std::make_shared<Virus<TSeq>>(*p.virus);
        virus->set_agent(this);
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
    
}

template<typename TSeq>
inline Agent<TSeq> & Agent<TSeq>::operator=(
    const Agent<TSeq> & other_agent
) 
{

    model = other_agent.model;

    n_neighbors = other_agent.n_neighbors;
    if (neighbors != nullptr)
    {
        delete neighbors;
        delete neighbors_locations;
    }

    if (other_agent.n_neighbors > 0u)
    {
        neighbors = new std::vector< size_t >(other_agent.n_neighbors);
        neighbors_locations = new std::vector< size_t >(other_agent.n_neighbors);
    }
    
    entities = other_agent.entities;
    entities_locations = other_agent.entities_locations;
    n_entities = other_agent.n_entities;

    state              = other_agent.state;
    state_prev         = other_agent.state_prev;
    state_last_changed = other_agent.state_last_changed;
    id                  = other_agent.id;
    
    if (other_agent.virus != nullptr)
    {
        virus = std::make_shared<Virus<TSeq>>(*other_agent.virus);
        virus->set_agent(this);
    } else
        virus = nullptr;
    
    n_tools             = other_agent.n_tools;
    for (size_t i = 0u; i < n_tools; ++i)
    {
        tools[i] = std::make_shared<Tool<TSeq>>(*other_agent.tools[i]);
        tools[i]->set_agent(this, i);
    }
    
    return *this;
    
}

template<typename TSeq>
inline Agent<TSeq>::~Agent()
{

    if (neighbors != nullptr)
    {
        delete neighbors;
        delete neighbors_locations;
    }

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

    CHECK_COALESCE_(state_new, tool->state_init, state);
    CHECK_COALESCE_(queue, tool->queue_init, Queue<TSeq>::NoOne);

    model->events_add(
        this, nullptr, tool, nullptr, state_new, queue, default_add_tool<TSeq>, -1, -1
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
inline void Agent<TSeq>::set_virus(
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

    CHECK_COALESCE_(state_new, virus->state_init, state);
    CHECK_COALESCE_(queue, virus->queue_init, Queue<TSeq>::NoOne);

    model->events_add(
        this, virus, nullptr, nullptr, state_new, queue, default_add_virus<TSeq>, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::set_virus(
    Virus<TSeq> virus,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    VirusPtr<TSeq> virus_ptr = std::make_shared< Virus<TSeq> >(virus);
    set_virus(virus_ptr, model, state_new, queue);
}

template<typename TSeq>
inline void Agent<TSeq>::add_entity(
    Entity<TSeq> & entity,
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    CHECK_COALESCE_(state_new, entity.state_init, state);
    CHECK_COALESCE_(queue, entity.queue_init, Queue<TSeq>::NoOne);

    if (model != nullptr)
    {

        model->events_add(
            this, nullptr, nullptr, &entity, state_new, queue, default_add_entity<TSeq>, -1, -1
        );

    }
    else // If no model is passed, then we assume that we only need to add the
         // model entity
    {

        Event<TSeq> a(
                this, nullptr, nullptr, &entity, state_new, queue, default_add_entity<TSeq>,
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

    CHECK_COALESCE_(state_new, tools[tool_idx]->state_post, state);
    CHECK_COALESCE_(queue, tools[tool_idx]->queue_post, Queue<TSeq>::NoOne);

    if (tool_idx >= n_tools)
        throw std::range_error(
            "The Tool you want to remove is out of range. This Agent only has " +
            std::to_string(n_tools) + " tools."
        );

    model->events_add(
        this, nullptr, tools[tool_idx], nullptr, state_new, queue, default_rm_tool<TSeq>, -1, -1
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

    model->events_add(
        this, nullptr, tool, nullptr, state_new, queue, default_rm_tool<TSeq>, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (virus == nullptr)
        throw std::logic_error(
            "There is no virus to remove here!"
        );

    CHECK_COALESCE_(state_new, virus->state_post, state);
    CHECK_COALESCE_(queue, virus->queue_post, Queue<TSeq>::Everyone);

    model->events_add(
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

    CHECK_COALESCE_(state_new, model->get_entity(entity_idx).state_post, state);
    CHECK_COALESCE_(queue, model->get_entity(entity_idx).queue_post, Queue<TSeq>::NoOne);

    model->events_add(
        this,
        nullptr,
        nullptr,
        &model->get_entity(entity_idx),
        state_new,
        queue, 
        default_rm_entity<TSeq>,
        entities_locations[entity_idx],
        entity_idx
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
        if (static_cast<int>(entities[i]) == entity.get_id())
        {
            entity_idx = i;
            break;
        }
    }

    if (entity_idx == -1)
        throw std::logic_error(
            std::string("The agent ") +
            std::to_string(id) +
            std::string(" is not associated with entity \"") +
            entity.get_name() +
            std::string("\".")
            );

    CHECK_COALESCE_(state_new, entity.state_post, state);
    CHECK_COALESCE_(queue, entity.queue_post, Queue<TSeq>::NoOne);

    model->events_add(
        this,
        nullptr,
        nullptr,
        &model->entities[entity.get_id()],
        state_new,
        queue, 
        default_rm_entity<TSeq>,
        entities_locations[entity_idx],
        entity_idx
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_agent_by_virus(
    Model<TSeq> * model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    CHECK_COALESCE_(state_new, virus->state_removed, state);
    CHECK_COALESCE_(queue, virus->queue_removed, Queue<TSeq>::Everyone);

    model->events_add(
        this, virus, nullptr, nullptr, state_new, queue,
        default_rm_virus<TSeq>, -1, -1
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
inline VirusPtr<TSeq> & Agent<TSeq>::get_virus() {
    return virus;
}

template<typename TSeq>
inline const VirusPtr<TSeq> & Agent<TSeq>::get_virus() const {
    return virus;
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

    virus->mutate();

}

template<typename TSeq>
inline void Agent<TSeq>::add_neighbor(
    Agent<TSeq> & p,
    bool check_source,
    bool check_target
) {
    // Can we find the neighbor?
    bool found = false;

    if (neighbors == nullptr)
    {
        neighbors = new std::vector< size_t >();
        neighbors_locations = new std::vector< size_t >();
    }

    if (check_source && neighbors)
    {

        for (auto & n: *neighbors)    
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

        neighbors_locations->push_back(p.get_n_neighbors());
        neighbors->push_back(p.get_id());
        n_neighbors++;

    }


    found = false;
    if (check_target && p.neighbors)
    {
       
        for (auto & n: *p.neighbors)
            if (static_cast<int>(n) == id)
            {
                found = true;
                break;
            }
    
    }

    if (!found)
    {

        if (p.neighbors == nullptr)
        {
            p.neighbors = new std::vector< size_t >();
            p.neighbors_locations = new std::vector< size_t >();
        }

        p.neighbors_locations->push_back(n_neighbors - 1);
        p.neighbors->push_back(id);
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

    if (n_this >= n_neighbors)
        throw std::range_error(
            "The neighbor you want to swap is out of range. This Agent only has " +
            std::to_string(n_neighbors) + " neighbors."
        );
    if (n_other >= other.n_neighbors)
        throw std::range_error(
            "The neighbor you want to swap is out of range. This Agent only has " +
            std::to_string(other.n_neighbors) + " neighbors."
        );

    // Getting the agents
    auto & pop = model->population;
    auto & neigh_this  = pop[(*neighbors)[n_this]];
    auto & neigh_other = pop[(*other.neighbors)[n_other]];

    // Getting the locations in the neighbors
    size_t loc_this_in_neigh = (*neighbors_locations)[n_this];
    size_t loc_other_in_neigh = (*other.neighbors_locations)[n_other];

    // Changing ids
    std::swap((*neighbors)[n_this], (*other.neighbors)[n_other]);

    if (!model->directed)
    {
        std::swap(
            (*neigh_this.neighbors)[loc_this_in_neigh],
            (*neigh_other.neighbors)[loc_other_in_neigh]
            );

        // Changing the locations
        std::swap((*neighbors_locations)[n_this], (*other.neighbors_locations)[n_other]);
        
        std::swap(
            (*neigh_this.neighbors_locations)[loc_this_in_neigh],
            (*neigh_other.neighbors_locations)[loc_other_in_neigh]
            );
    }

}

template<typename TSeq>
inline std::vector< Agent<TSeq> *> Agent<TSeq>::get_neighbors()
{
    std::vector< Agent<TSeq> * > res(n_neighbors, nullptr);
    for (size_t i = 0u; i < n_neighbors; ++i)
        res[i] = &model->population[(*neighbors)[i]];

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

    model->events_add(
        this, nullptr, nullptr, nullptr, new_state, queue,
        default_change_state<TSeq>, -1, -1
    );
    
    return;

}

template<typename TSeq>
inline const unsigned int & Agent<TSeq>::get_state() const {
    return state;
}

template<typename TSeq>
inline void Agent<TSeq>::reset()
{

    this->virus = nullptr;

    this->tools.clear();
    n_tools = 0u;

    this->entities.clear();
    this->entities_locations.clear();
    this->n_entities = 0u;

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
    if (virus->get_id() == static_cast<int>(t))
        return true;

    return false;
}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(std::string name) const
{
    
    if (virus->get_name() == name)
        return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(const Virus<TSeq> & virus) const
{

    return has_virus(virus.get_id());

}

template<typename TSeq>
inline bool Agent<TSeq>::has_entity(epiworld_fast_uint t) const
{

    for (auto & entity : entities)
        if (entity == t)
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_entity(std::string name) const
{

    for (auto & entity : entities)
        if (model->get_entity(entity).get_name() == name)
            return true;

    return false;

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
            "Agent: %i, state: %s (%i), Has virus: %s, NTools: %ii NNeigh: %i\n",
            static_cast<int>(id),
            model->states_labels[state].c_str(),
            static_cast<int>(state),
            virus == nullptr ? std::string("no").c_str() : std::string("yes").c_str(),
            static_cast<int>(n_tools),
            static_cast<int>(n_neighbors)
        );
    }
    else {

        printf_epiworld("Information about agent id %i\n",
            static_cast<int>(this->id));
        printf_epiworld("  State        : %s (%i)\n",
            model->states_labels[state].c_str(), static_cast<int>(state));
        printf_epiworld("  Has virus    : %s\n", virus == nullptr ?
            std::string("no").c_str() : std::string("yes").c_str());
        printf_epiworld("  Tool count   : %i\n", static_cast<int>(n_tools));
        printf_epiworld("  Neigh. count : %i\n", static_cast<int>(n_neighbors));

        size_t nfeats = model->get_agents_data_ncols();
        if (nfeats > 0)
        {

            printf_epiworld(
                "This model includes features (%i): [ ",
                static_cast<int>(nfeats)
                );

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
    if (n_entities == 0)
        throw std::range_error("Agent id " + std::to_string(id) + " has no entities.");

    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->get_entity(entities[i]);
}

template<typename TSeq>
inline Entity<TSeq> & Agent<TSeq>::get_entity(size_t i)
{
    if (n_entities == 0)
        throw std::range_error("Agent id " + std::to_string(id) + " has no entities.");

    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->get_entity(entities[i]);
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
        ((virus == nullptr) && (other.virus != nullptr)) ||
            ((virus != nullptr) && (other.virus == nullptr)),
        "Agent:: virus don't match"
    )

    if ((virus != nullptr) && (other.virus != nullptr))
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *virus != *other.virus,
            "Agent:: virus doesn't match"
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