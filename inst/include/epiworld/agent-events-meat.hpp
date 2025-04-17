#ifndef EPIWORLD_AGENT_EVENTS_MEAT_HPP
#define EPIWORLD_AGENT_EVENTS_MEAT_HPP

template<typename TSeq>
inline void default_add_virus(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> *  p = a.agent;
    VirusPtr<TSeq> v = a.virus;
    
    m->get_db().record_transmission(
        v->get_agent() ? v->get_agent()->get_id() : -1,
        p->get_id(),
        v->get_id(),
        v->get_date() 
    );
    
    p->virus = std::make_shared< Virus<TSeq> >(*v);
    p->virus->set_date(m->today());
    p->virus->set_agent(p);

    // Change of state needs to be recorded and updated on the
    // tools.
    if (p->state_prev != p->state)
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, p->state);

        for (size_t i = 0u; i < p->n_tools; ++i)
            db.update_tool(p->tools[i]->get_id(), p->state_prev, p->state);
    }

    // Lastly, we increase the daily count of the virus
    #ifdef EPI_DEBUG
    m->get_db().today_virus.at(v->get_id()).at(p->state)++;
    #else
    m->get_db().today_virus[v->get_id()][p->state]++;
    #endif

}

template<typename TSeq>
inline void default_add_tool(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> t = a.tool;
    
    // Update tool accounting
    p->n_tools++;
    size_t n_tools = p->n_tools;

    if (n_tools <= p->tools.size())
        p->tools[n_tools - 1] = std::make_shared< Tool<TSeq> >(*t);
    else
        p->tools.push_back(std::make_shared< Tool<TSeq> >(*t));

    n_tools--;

    p->tools[n_tools]->set_date(m->today());
    p->tools[n_tools]->set_agent(p, n_tools);

    // Change of state needs to be recorded and updated on the
    // tools.
    if (p->state_prev != p->state)
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, p->state);

        if (p->virus)
            db.update_virus(p->virus->get_id(), p->state_prev, p->state);
    }

    m->get_db().today_tool[t->get_id()][p->state]++;


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
    if (p->state_prev != p->state)
    {
        auto & db = model->get_db();
        db.update_state(p->state_prev, p->state);

        for (size_t i = 0u; i < p->n_tools; ++i)
            db.update_tool(p->tools[i]->get_id(), p->state_prev, p->state);
    }

    // The counters of the virus only needs to decrease
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

    if (--p->n_tools > 0)
    {
        p->tools[p->n_tools]->pos_in_agent = t->pos_in_agent;
        std::swap(
            p->tools[t->pos_in_agent],
            p->tools[p->n_tools]
            );
    }

    // Change of state needs to be recorded and updated on the
    // tools.
    if (p->state_prev != p->state)
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, p->state);

        if (p->virus)
            db.update_virus(p->virus->get_id(), p->state_prev, p->state);
    }

    // Lastly, we increase the daily count of the tool
    #ifdef EPI_DEBUG
    m->get_db().today_tool.at(t->get_id()).at(p->state_prev)--;
    #else
    m->get_db().today_tool[t->get_id()][p->state_prev]--;
    #endif

    return;

}

template<typename TSeq>
inline void default_change_state(Event<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;

    if (p->state_prev != p->state)
    {
        auto & db = m->get_db();
        db.update_state(p->state_prev, p->state);

        if (p->virus)
            db.update_virus(p->virus->get_id(), p->state_prev, p->state);

        for (size_t i = 0u; i < p->n_tools; ++i)
            db.update_tool(p->tools[i]->get_id(), p->state_prev, p->state);

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
    if (++p->n_entities <= p->entities.size())
    {

        p->entities[p->n_entities - 1]           = e->get_id();
        p->entities_locations[p->n_entities - 1] = e->n_agents;

    } else
    {
        p->entities.push_back(e->get_id());
        p->entities_locations.push_back(e->n_agents);
    }

    // Adding the agent to the entity
    // Adding the entity to the agent
    if (++e->n_agents <= e->agents.size())
    {

        e->agents[e->n_agents - 1]          = p->get_id();
        // Adjusted by '-1' since the list of entities in the agent just grew.
        e->agents_location[e->n_agents - 1] = p->n_entities - 1;

    } else
    {
        e->agents.push_back(p->get_id());
        e->agents_location.push_back(p->n_entities - 1);
    }

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

    if (--p->n_entities > 0)
    {

        // When we move the end entity to the new location, the 
        // moved entity needs to reflect the change, i.e., where the
        // entity will now be located in the agent
        size_t agent_location_in_last_entity  =
            p->entities_locations[p->n_entities];

        Entity<TSeq> * last_entity =
            &m->get_entity(p->entities[p->n_entities]); ///< Last entity of the agent

        // The end entity will be located where the removed was
        last_entity->agents_location[agent_location_in_last_entity] =
            idx_entity_in_agent;

        // We now make the swap
        std::swap(
            p->entities[p->n_entities],
            p->entities[idx_entity_in_agent]
        );

    }

    if (--e->n_agents > 0)
    {

        // When we move the end agent to the new location, the 
        // moved agent needs to reflect the change, i.e., where the
        // agent will now be located in the entity
        size_t entity_location_in_last_agent = e->agents_location[e->n_agents];
        
        Agent<TSeq> * last_agent  =
            &m->get_agents()[e->agents[e->n_agents]]; ///< Last agent of the entity

        // The end entity will be located where the removed was
        last_agent->entities_locations[entity_location_in_last_agent] =
            idx_agent_in_entity;

        // We now make the swap
        std::swap(
            e->agents[e->n_agents],
            e->agents[idx_agent_in_entity]
        );

    }

    // Setting the date of the last removal
    // e->date_last_add_or_remove = m->today();

    return;

}
#endif