#ifndef EPIWORLD_AGENT_ACTIONS_MEAT_HPP
#define EPIWORLD_AGENT_ACTIONS_MEAT_HPP

template<typename TSeq>
inline void default_add_virus(Action<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> *  p = a.agent;
    VirusPtr<TSeq> v = a.virus;

    CHECK_COALESCE_(a.new_state, v->state_init, p->get_state())
    CHECK_COALESCE_(a.queue, v->queue_init, 1)

    // Has a agent? If so, we need to register the transmission
    if (v->get_agent())
    {

        // ... only if not the same agent
        if (v->get_agent()->get_id() != p->get_id())
            m->get_db().record_transmission(
                v->get_agent()->get_id(),
                p->get_id(),
                v->get_id(),
                v->get_date() 
            );

    }
    
    // Update virus accounting
    p->n_viruses++;
    size_t n_viruses = p->n_viruses;

    if (n_viruses <= p->viruses.size())
        p->viruses[n_viruses - 1] = std::make_shared< Virus<TSeq> >(*v);
    else
        p->viruses.push_back(std::make_shared< Virus<TSeq> >(*v));

    n_viruses--;

    // Notice that both agent and date can be changed in this case
    // as only the sequence is a shared_ptr itself.
    #ifdef EPI_DEBUG
    if (n_viruses >= p->viruses.size())
    {
        throw std::logic_error(
            "[epi-debug]::default_add_virus Index for new virus out of range."
            );
    }
    #endif
    p->viruses[n_viruses]->set_agent(p, n_viruses);
    p->viruses[n_viruses]->set_date(m->today());

    #ifdef EPI_DEBUG
    m->get_db().today_virus.at(v->get_id()).at(p->state)++;
    #else
    m->get_db().today_virus[v->get_id()][p->state]++;
    #endif

}

template<typename TSeq>
inline void default_add_tool(Action<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> t = a.tool;

    CHECK_COALESCE_(a.new_state, t->state_init, p->get_state())
    CHECK_COALESCE_(a.queue, t->queue_init, Queue<TSeq>::NoOne)
    
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

    m->get_db().today_tool[t->get_id()][p->state]++;

}

template<typename TSeq>
inline void default_rm_virus(Action<TSeq> & a, Model<TSeq> * model)
{

    Agent<TSeq> * p    = a.agent;    
    VirusPtr<TSeq> & v = a.agent->viruses[a.virus->pos_in_agent];
    
    CHECK_COALESCE_(a.new_state, v->state_post, p->get_state())
    CHECK_COALESCE_(a.queue, v->queue_post, -Queue<TSeq>::Everyone)

    if (--p->n_viruses > 0)
    {
        // The new virus will change positions
        p->viruses[p->n_viruses]->pos_in_agent = v->pos_in_agent;
        std::swap(
            p->viruses[p->n_viruses],   // Moving to the end
            p->viruses[v->pos_in_agent] // Moving to the beginning
            );
    }
    
    // Calling the virus action over the removed virus
    v->post_recovery(model);

    return;

}

template<typename TSeq>
inline void default_rm_tool(Action<TSeq> & a, Model<TSeq> * /*m*/)
{

    Agent<TSeq> * p   = a.agent;    
    ToolPtr<TSeq> & t = a.agent->tools[a.tool->pos_in_agent];

    CHECK_COALESCE_(a.new_state, t->state_post, p->get_state())
    CHECK_COALESCE_(a.queue, t->queue_post, Queue<TSeq>::NoOne)

    if (--p->n_tools > 0)
    {
        p->tools[p->n_tools]->pos_in_agent = t->pos_in_agent;
        std::swap(
            p->tools[t->pos_in_agent],
            p->tools[p->n_tools]
            );
    }

    return;

}

template<typename TSeq>
inline void default_add_entity(Action<TSeq> & a, Model<TSeq> *)
{

    Agent<TSeq> *  p = a.agent;
    Entity<TSeq> * e = a.entity;

    CHECK_COALESCE_(a.new_state, e->state_post, p->get_state())
    CHECK_COALESCE_(a.queue, e->queue_post, Queue<TSeq>::NoOne)

    // Checking the agent and the entity are not linked
    if ((p->get_n_entities() > 0) && (e->size() > 0))
    {

        if (p->get_n_entities() > e->size()) // Slower search through the agent
        {
            for (size_t i = 0u; i < e->size(); ++i)
                if(e->operator[](i)->get_id() == p->get_id())
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
inline void default_rm_entity(Action<TSeq> & a, Model<TSeq> * m)
{
    
    Agent<TSeq> *  p = a.agent;    
    Entity<TSeq> * e = a.entity;
    size_t idx_agent_in_entity = a.idx_agent;
    size_t idx_entity_in_agent = a.idx_object;

    CHECK_COALESCE_(a.new_state, e->state_post, p->get_state())
    CHECK_COALESCE_(a.queue, e->queue_post, Queue<TSeq>::NoOne)

    if (--p->n_entities > 0)
    {

        // When we move the end entity to the new location, the 
        // moved entity needs to reflect the change, i.e., where the
        // entity will now be located in the agent
        size_t agent_location_in_last_entity  = p->entities_locations[p->n_entities];
        Entity<TSeq> * last_entity = &m->get_entities()[p->entities[p->n_entities]]; ///< Last entity of the agent

        // The end entity will be located where the removed was
        last_entity->agents_location[agent_location_in_last_entity] = idx_entity_in_agent;

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
        Agent<TSeq> * last_agent  = &m->get_agents()[e->agents[e->n_agents]]; ///< Last agent of the entity

        // The end entity will be located where the removed was
        last_agent->entities_locations[entity_location_in_last_agent] = idx_agent_in_entity;

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