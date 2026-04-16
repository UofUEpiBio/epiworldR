#ifndef EPIWORLD_PERSON_MEAT_HPP
#define EPIWORLD_PERSON_MEAT_HPP

#include <vector>
#include <string>
#include "config.hpp"
#include "epiworld-macros.hpp"
#include "agent-bones.hpp"
#include "agent-events-meat.hpp"

template<typename TSeq>
inline Agent<TSeq>::Agent(Agent<TSeq> && p) :
    neighbors(std::move(p.neighbors)),
    neighbors_locations(std::move(p.neighbors_locations)),
    n_neighbors(p.n_neighbors),
    entities(std::move(p.entities)),
    state(p.state),
    state_prev(p.state_prev), 
    state_last_changed(p.state_last_changed),
    id(p.id),
    tools(std::move(p.tools)) /// Needs to be adjusted
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

// Copy constructor
template<typename TSeq>
inline Agent<TSeq>::Agent(const Agent<TSeq> & p) :
    neighbors(nullptr),
    neighbors_locations(nullptr),
    n_neighbors(p.n_neighbors),
    entities(p.entities)
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
        virus = std::shared_ptr<Virus<TSeq>>(p.virus->clone_ptr());
        virus->set_agent(this);
    }
    

    tools.clear();
    tools.reserve(p.get_n_tools());
    for (size_t i = 0u; i < p.get_n_tools(); ++i)
    {
        tools.emplace_back(std::shared_ptr<Tool<TSeq>>(p.tools[i]->clone_ptr()));
        tools.back()->set_agent(this, i);

    }
    
}

template<typename TSeq>
inline Agent<TSeq> & Agent<TSeq>::operator=(
    const Agent<TSeq> & other_agent
)
{

    n_neighbors = other_agent.n_neighbors;
    if (neighbors != nullptr)
    {
        delete neighbors;
        delete neighbors_locations;
    }

    if (other_agent.n_neighbors > 0u)
    {
        neighbors = new std::vector< size_t >(*other_agent.neighbors);
        neighbors_locations = new std::vector< size_t >(*other_agent.neighbors_locations);
    }
    else 
    {
        neighbors = nullptr;
        neighbors_locations = nullptr;
    }
    
    entities = other_agent.entities;

    state              = other_agent.state;
    state_prev         = other_agent.state_prev;
    state_last_changed = other_agent.state_last_changed;
    id                  = other_agent.id;
    
    if (other_agent.virus != nullptr)
    {
        virus = std::shared_ptr<Virus<TSeq>>(other_agent.virus->clone_ptr());
        virus->set_agent(this);
    } else
        virus = nullptr;
    
    
    tools.clear();
    tools.reserve(other_agent.get_n_tools());
    for (size_t i = 0u; i < other_agent.get_n_tools(); ++i)
    {
        tools.emplace_back(std::shared_ptr<Tool<TSeq>>(other_agent.tools[i]->clone_ptr()));
        tools.back()->set_agent(this, i);
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
    Model<TSeq> & model,
    const Tool<TSeq> & tool,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
) {

    // Checking the tool exists
    if (tool.get_id() >= static_cast<int>(model.get_db().get_n_tools()))
        throw std::range_error("The tool with id: " + std::to_string(tool.get_id()) +
            " has not been registered. There are only " + std::to_string(model.get_n_tools()) +
            " included in the model.");

    ToolPtr<TSeq> tool_ptr = std::shared_ptr<Tool<TSeq>>(tool.clone_ptr());

    model._add_event(
        this, nullptr, tool_ptr, nullptr, state_new, queue, EventAction::AddTool
        );

}

template<typename TSeq>
inline void Agent<TSeq>::set_virus(
    Model<TSeq> & model,
    const Virus<TSeq> & virus,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    
        // Checking the virus exists
    if (virus.get_id() >= static_cast<int>(model.get_db().get_n_viruses()))
        throw std::range_error("The virus with id: " + std::to_string(virus.get_id()) +
            " has not been registered. There are only " + std::to_string(model.get_n_viruses()) +
            " included in the model.");

    if (state_new == -99)
        virus.get_state(&state_new, nullptr, nullptr);

    if (queue == -99)
        virus.get_queue(&queue, nullptr, nullptr);

    VirusPtr<TSeq> virus_ptr = std::shared_ptr<Virus<TSeq>>(virus.clone_ptr());

    model._add_event(
        this, virus_ptr, nullptr, nullptr, state_new, queue, EventAction::AddVirus
        );

}

template<typename TSeq>
inline void Agent<TSeq>::add_entity(
    Model<TSeq> & model,
    Entity<TSeq> & entity,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    model._add_event(
        this, nullptr, nullptr, &entity, state_new, queue, EventAction::AddEntity
    );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    Model<TSeq> & model,
    epiworld_fast_uint tool_idx,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool_idx >= tools.size())
        throw std::range_error(
            "The Tool you want to remove is out of range. This Agent only has " +
            std::to_string(tools.size()) + " tools."
        );

    model._add_event(
        this, nullptr, tools[tool_idx], nullptr, state_new, queue, EventAction::RemoveTool
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    Model<TSeq> & model,
    ToolPtr<TSeq> & tool,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model._add_event(
        this, nullptr, tool, nullptr, state_new, queue, EventAction::RemoveTool
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    Model<TSeq> & model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (virus == nullptr)
        throw std::logic_error(
            "There is no virus to remove here!"
        );

    if (state_new == -99)
        virus->get_state(nullptr, &state_new, nullptr);

    if (queue == -99)
        virus->get_queue(nullptr, &queue, nullptr);

    model._add_event(
        this, virus, nullptr, nullptr,
        state_new,
        queue,
        EventAction::RemoveVirus
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    Model<TSeq> & model,
    epiworld_fast_uint entity_idx,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (entity_idx >= entities.size())
        throw std::range_error(
            "The Entity you want to remove is out of range. This Agent only has " +
            std::to_string(entities.size()) + " entity(ies)."
        );
    else if (entities.size() == 0u)
        throw std::logic_error(
            "There is no entity to remove here!"
        );

    model._add_event(
        this,
        nullptr,
        nullptr,
        &model.get_entity(entities[entity_idx]),
        state_new,
        queue,
        EventAction::RemoveEntity
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    Model<TSeq> & model,
    Entity<TSeq> & entity,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    // Looking for entity location in the agent
    bool found = std::find(
        entities.begin(), entities.end(),
        static_cast<size_t>(entity.get_id())
    ) != entities.end();

    if (!found)
        throw std::logic_error(
            std::string("The agent ") +
            std::to_string(id) +
            std::string(" is not associated with entity \"") +
            entity.get_name() +
            std::string("\".")
            );

    model._add_event(
        this,
        nullptr,
        nullptr,
        &model.get_entity(entity.get_id()),
        state_new,
        queue,
        EventAction::RemoveEntity
    );
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> & model
) {

    return model.susceptibility_reduction_mixer(this, v);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> & model
) {
    return model.transmission_reduction_mixer(this, v);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> & v,
    Model<TSeq> & model
) {
    return model.recovery_enhancer_mixer(this, v);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_death_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> & model
) {
    return model.death_reduction_mixer(this, v);
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
inline std::vector< ToolPtr<TSeq> > Agent<TSeq>::get_tools() {
    return tools;
}

template<typename TSeq>
inline const std::vector< ToolPtr<TSeq> > Agent<TSeq>::get_tools() const {
    return tools;
}

template<typename TSeq>
inline ToolPtr<TSeq> & Agent<TSeq>::get_tool(int i)
{
    return tools.at(i);
}

template<typename TSeq>
inline ToolPtr<TSeq> & Agent<TSeq>::get_tool(std::string name)
{
    for (auto & tool : tools)
        if (tool->get_name() == name)
            return tool;

    throw std::logic_error(
        "The agent does not have a tool with name: " + name
    );

}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_tools() const noexcept
{
    return tools.size();
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
    size_t n_other,
    Model<TSeq> & model
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
    auto & pop = model.population;
    auto & neigh_this  = pop[(*neighbors)[n_this]];
    auto & neigh_other = pop[(*other.neighbors)[n_other]];

    // Getting the locations in the neighbors
    size_t loc_this_in_neigh = (*neighbors_locations)[n_this];
    size_t loc_other_in_neigh = (*other.neighbors_locations)[n_other];

    // Changing ids
    std::swap((*neighbors)[n_this], (*other.neighbors)[n_other]);

    if (!model.directed)
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
inline std::vector< Agent<TSeq> *> Agent<TSeq>::get_neighbors(Model<TSeq> & model)
{
    std::vector< Agent<TSeq> * > res(n_neighbors, nullptr);
    for (size_t i = 0u; i < n_neighbors; ++i)
        res[i] = &model.population[(*neighbors)[i]];
    return res;
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_neighbors() const
{
    return n_neighbors;
}

template<typename TSeq>
inline void Agent<TSeq>::change_state(
    Model<TSeq> & model,
    epiworld_fast_uint new_state,
    epiworld_fast_int queue
    )
{

    model._add_event(
        this, nullptr, nullptr, nullptr, new_state, queue,
        EventAction::ChangeState
    );

    return;

}

template<typename TSeq>
inline unsigned int Agent<TSeq>::get_state() const {
    return state;
}

template<typename TSeq>
inline unsigned int Agent<TSeq>::get_state_prev() const {
    return state_prev;
}

template<typename TSeq>
inline int Agent<TSeq>::get_state_last_changed() const {
    return state_last_changed;
}

template<typename TSeq>
inline void Agent<TSeq>::reset()
{

    this->virus = nullptr;

    this->tools.clear();
    decltype(this->tools)().swap(this->tools);

    this->entities.clear();
    decltype(this->entities)().swap(this->entities);

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
inline bool Agent<TSeq>::has_tool(std::string_view name) const
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
inline bool Agent<TSeq>::has_virus(std::string_view name) const
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

    for (auto & entity_id : entities)
        if (entity_id == static_cast<size_t>(t))
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_entity(std::string_view name, const Model<TSeq> & model) const
{

    for (auto & entity_id : entities)
        if (model.get_entity(entity_id).get_name() == name)
            return true;

    return false;

}

template<typename TSeq>
inline void Agent<TSeq>::print(
    Model<TSeq> & model,
    bool compressed
    ) const
{

    if (compressed)
    {
        printf_epiworld(
            "Agent: %i, state: %s (%i), Has virus: %s, NTools: %ii NNeigh: %i\n",
            static_cast<int>(id),
            model.states_labels[state].c_str(),
            static_cast<int>(state),
            virus == nullptr ? std::string("no").c_str() : std::string("yes").c_str(),
            static_cast<int>(tools.size()),
            static_cast<int>(n_neighbors)
        );
    }
    else {

        printf_epiworld("Information about agent id %i\n",
            static_cast<int>(this->id));
        printf_epiworld("  State        : %s (%i)\n",
            model.states_labels[state].c_str(), static_cast<int>(state));
        printf_epiworld("  Has virus    : %s\n", virus == nullptr ?
            std::string("no").c_str() : std::string("yes").c_str());
        printf_epiworld("  Tool count   : %i\n", static_cast<int>(tools.size()));
        printf_epiworld("  Neigh. count : %i\n", static_cast<int>(n_neighbors));

        size_t nfeats = model.get_agents_data_ncols();
        if (nfeats > 0)
        {

            printf_epiworld(
                "This model includes features (%i): [ ",
                static_cast<int>(nfeats)
                );

            int max_to_show = static_cast<int>((nfeats > 10)? 10 : nfeats);

            for (int k = 0; k < max_to_show; ++k)
            {
                printf_epiworld("%.2f", this->operator()(k, model));

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
inline double & Agent<TSeq>::operator()(size_t j, Model<TSeq> & model)
{

    if (model.agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model.agents_data + j * model.size() + id);

}

template<typename TSeq>
inline double Agent<TSeq>::operator()(size_t j, const Model<TSeq> & model) const
{

    if (model.agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model.agents_data + j * model.size() + id);

}

template<typename TSeq>
inline const std::vector< size_t > & Agent<TSeq>::get_entities() const
{
    return entities;
}

template<typename TSeq>
inline const Entity<TSeq> & Agent<TSeq>::get_entity(size_t i, const Model<TSeq> & model) const
{
    if (entities.size() == 0)
        throw std::range_error("Agent id " + std::to_string(id) + " has no entities.");

    if (i >= entities.size())
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model.get_entity(
        #ifdef EPI_DEBUG
        entities.at(i)
        #else
        entities[i]
        #endif
    );
}

template<typename TSeq>
inline Entity<TSeq> & Agent<TSeq>::get_entity(size_t i, Model<TSeq> & model)
{
    if (entities.size() == 0)
        throw std::range_error("Agent id " + std::to_string(id) + " has no entities.");

    if (i >= entities.size())
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model.get_entity(
        #ifdef EPI_DEBUG
        entities.at(i)
        #else
        entities[i]
        #endif
    );
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_entities() const
{
    return entities.size();
}

template<typename TSeq>
inline bool Agent<TSeq>::operator==(const Agent<TSeq> & other) const
{

    // Checking the address
    if (this == &other)
        return true;

    EPI_DEBUG_FAIL_AT_TRUE(
        n_neighbors != other.n_neighbors,
        "Agent:: n_eighbors don't match"
        )

    
    for (size_t i = 0u; i < n_neighbors; ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            (*neighbors)[i] != (*other.neighbors)[i],
            "Agent:: neighbor[i] don't match"
        )
    }
    
    EPI_DEBUG_FAIL_AT_TRUE(
        entities.size() != other.entities.size(),
        "Agent:: n_entities don't match"
        )


    for (size_t i = 0u; i < entities.size(); ++i)
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
    
    EPI_DEBUG_FAIL_AT_TRUE(
        tools.size() != other.tools.size(),
        "Agent:: n_tools don't match"
    )

    for (size_t i = 0u; i < tools.size(); ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            *tools[i] != *other.tools[i],
            "Agent:: tools[i] don't match"
        )
         
    }   
    
    return true;
    
}

#endif
