#ifndef EPIWORLD_AGENT_EVENTS_MEAT_HPP
#define EPIWORLD_AGENT_EVENTS_MEAT_HPP


template<typename TSeq>
inline void default_add_virus(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> *  p = a.agent;
    VirusPtr<TSeq> & v = a.virus;
    
    m->get_db().record_transmission(
        v->get_agent() ? v->get_agent()->get_id() : -1,
        p->get_id(),
        v->get_id(),
        v->get_date() 
    );
    
    p->virus = std::move(v);
    p->virus->set_date(m->today());
    p->virus->set_agent(p);

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, a.new_state);

        // For tool counts, use current state (p->state) not state_prev
        // because state_prev may be stale if multiple changes occurred today
        for (size_t i = 0u; i < p->tools.size(); ++i)
            db.update_tool(
                p->tools[i]->get_id(),
                p->state,
                a.new_state
            );
    }

    // Lastly, we increase the daily count of the virus
    #ifdef EPI_DEBUG
    m->get_db().today_virus.at(p->virus->get_id()).at(
        a.new_state != -99 ? a.new_state : p->state
    )++;
    #else
    m->get_db().today_virus[p->virus->get_id()][
        a.new_state != -99 ? a.new_state : p->state
    ]++;
    #endif

}

template<typename TSeq>
inline void default_add_tool(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> & t = a.tool;
    
    // Update tool accounting
    p->tools.emplace_back(std::move(t));

    size_t last_tool_idx = p->tools.size() - 1u;

    p->tools[last_tool_idx]->set_date(m->today());
    p->tools[last_tool_idx]->set_agent(p, last_tool_idx);

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && static_cast<int>(p->state) != a.new_state)
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, a.new_state);

        // For virus counts, use current state (p->state) not state_prev
        // because state_prev may be stale if multiple changes occurred today
        if (p->virus)
            db.update_virus(
                p->virus->get_id(),
                p->state,
                a.new_state
            );
    }

    m->get_db().today_tool[p->tools.back()->get_id()][
        a.new_state != -99 ? a.new_state : p->state
    ]++;


}

template<typename TSeq>
inline void default_rm_virus(Event<TSeq> & a, Model<TSeq> * model)
{

    Agent<TSeq> * p    = a.agent;
    VirusPtr<TSeq> & v = a.virus;

    // Calling the virus action over the removed virus
    v->post_recovery(model);

    p->virus = nullptr;

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
        auto & db = model->get_db();
        db.update_state(p->state_prev, a.new_state);

        // For tool counts, use current state (p->state) not state_prev
        // because state_prev may be stale if multiple changes occurred today
        for (size_t i = 0u; i < p->tools.size(); ++i)
            db.update_tool(
                p->tools[i]->get_id(),
                p->state,
                a.new_state
            );
    }

    // The counters of the virus only needs to decrease.
    // We use the previous state of the agent as that was
    // the state when the virus was added.
    #ifdef EPI_DEBUG
    model->get_db().today_virus.at(v->get_id()).at(p->state_prev)--;
    #else
    model->get_db().today_virus[v->get_id()][p->state_prev]--;
    #endif

    
    return;

}

template<typename TSeq>
inline void default_rm_tool(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p   = a.agent;    
    ToolPtr<TSeq> & t = a.agent->tools[a.tool->pos_in_agent];

    // Capture the tool ID before any swaps or removals
    int removed_tool_id = t->get_id();

    if (p->tools.size() > 1u)
    {
        p->tools.back()->pos_in_agent = t->pos_in_agent;
        std::swap(
            p->tools[t->pos_in_agent],
            p->tools.back()
            );
    }

    // Remove the last element (the one we want to remove)
    p->tools.pop_back();

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, a.new_state);

        // For virus counts, use current state (p->state) not state_prev
        // because state_prev may be stale if multiple changes occurred today
        if (p->virus)
            db.update_virus(
                p->virus->get_id(),
                p->state,
                a.new_state
            );
    }

    // Lastly, we increase the daily count of the tool.
    // Like rm_virus, we use the previous state of the agent
    // as that was the state when the tool was added.
    #ifdef EPI_DEBUG
    m->get_db().today_tool.at(removed_tool_id).at(p->state_prev)--;
    #else
    m->get_db().today_tool[removed_tool_id][p->state_prev]--;
    #endif

    return;

}

template<typename TSeq>
inline void default_change_state(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;

    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, a.new_state);

        // For virus and tool counts, use current state (p->state) not state_prev
        // because state_prev may be stale if multiple changes occurred today
        if (p->virus)
            db.update_virus(
                p->virus->get_id(), p->state, a.new_state
            );

        for (size_t i = 0u; i < p->tools.size(); ++i)
            db.update_tool(
                p->tools[i]->get_id(),
                p->state,
                a.new_state
            );

    }

}

template<typename TSeq>
inline void default_add_entity(Event<TSeq> & a, Model<TSeq> *)
{

    Agent<TSeq> *  p = a.agent;
    Entity<TSeq> * e = a.entity;

    // Checking the agent and the entity are not linked
    if ((p->get_n_entities() > 0) && (e->size() > 0))
    {

        if (p->get_n_entities() > e->size()) // Slower search through the agent
        {
            for (size_t i = 0u; i < e->size(); ++i)
                if(static_cast<int>(e->operator[](i)) == p->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }
        else                                 // Slower search through the entity
        {
            for (size_t i = 0u; i < p->get_n_entities(); ++i)
                if(p->get_entity(i).get_id() == e->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }

        // It means that agent and entity were not associated.
    }

    // Adding the entity to the agent
    p->entities.push_back(e->get_id());
    p->entities_locations.push_back(e->agents.size());

    // Adding the agent to the entity
    e->agents.push_back(p->get_id());
    e->agents_location.push_back(p->entities.size() - 1);

    // Today was the last modification
    // e->date_last_add_or_remove = m->today();
    
}

template<typename TSeq>
inline void default_rm_entity(Event<TSeq> & a, Model<TSeq> * m)
{
    
    Agent<TSeq> *  p = a.agent;    
    Entity<TSeq> * e = a.entity;
    size_t idx_agent_in_entity = a.idx_agent;
    size_t idx_entity_in_agent = a.idx_object;

    if (p->entities.size() > 1u)
    {

        // When we move the end entity to the new location, the 
        // moved entity needs to reflect the change, i.e., where the
        // entity will now be located in the agent
        size_t agent_location_in_last_entity  =
            p->entities_locations.back();

        Entity<TSeq> * last_entity =
            &m->get_entity(p->entities.back()); ///< Last entity of the agent

        // The end entity will be located where the removed was
        last_entity->agents_location[agent_location_in_last_entity] =
            idx_entity_in_agent;

        // We now make the swap
        std::swap(
            p->entities.back(),
            p->entities[idx_entity_in_agent]
        );

        std::swap(
            p->entities_locations.back(),
            p->entities_locations[idx_entity_in_agent]
        );

    }

    // Remove the last element (the one we want to remove)
    p->entities.pop_back();
    p->entities_locations.pop_back();

    if (e->agents.size() > 1u)
    {

        // When we move the end agent to the new location, the 
        // moved agent needs to reflect the change, i.e., where the
        // agent will now be located in the entity
        size_t entity_location_in_last_agent = e->agents_location.back();
        
        Agent<TSeq> * last_agent  =
            &m->get_agents()[e->agents.back()]; ///< Last agent of the entity

        // The end entity will be located where the removed was
        last_agent->entities_locations[entity_location_in_last_agent] =
            idx_agent_in_entity;

        // We now make the swap
        std::swap(
            e->agents.back(),
            e->agents[idx_agent_in_entity]
        );

        std::swap(
            e->agents_location.back(),
            e->agents_location[idx_agent_in_entity]
        );

    }

    // Remove the last element (the one we want to remove)
    e->agents.pop_back();
    e->agents_location.pop_back();

    // Setting the date of the last removal
    // e->date_last_add_or_remove = m->today();

    return;

};

#endif