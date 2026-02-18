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
        for (size_t i = 0u; i < p->n_tools; ++i)
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
    p->n_tools++;
    size_t n_tools = p->n_tools;

    if (n_tools <= p->tools.size())
        p->tools[n_tools - 1] = std::move(t);
    else
        p->tools.emplace_back(std::move(t));

    n_tools--;

    p->tools[n_tools]->set_date(m->today());
    p->tools[n_tools]->set_agent(p, n_tools);

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
        for (size_t i = 0u; i < p->n_tools; ++i)
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

        for (size_t i = 0u; i < p->n_tools; ++i)
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
            for (const Agent<TSeq> & agent: e->get_agents())
                if(agent.get_id() == p->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }
        else                                 // Slower search through the entity
        {
            for (const Entity<TSeq> & entity: p->get_entities())
                if(entity.get_id() == e->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }

        // It means that agent and entity were not associated.
    }

    // Adding the to agent and the entity
    p->entities.push_back(std::ref(*e));
    e->agents.push_back(std::ref(*p));
    
}

template<typename TSeq>
inline void default_rm_entity(Event<TSeq> & a, Model<TSeq> *)
{
    
    Agent<TSeq> &  p = *a.agent;    
    Entity<TSeq> & e = *a.entity;
    
    // Remove entity from agent's entity list
    p.entities.erase(
        std::remove_if(
            p.entities.begin(),
            p.entities.end(),
            [&e](const std::reference_wrapper<Entity<TSeq>> & entity_ref) {
                return entity_ref.get().get_id() == e.get_id();
            }
        ),
        p.entities.end()
    );

    // Remove agent from entity's agent list
    e.agents.erase(
        std::remove_if(
            e.agents.begin(),
            e.agents.end(),
            [&p](const std::reference_wrapper<Agent<TSeq>> & agent_ref) {
                return agent_ref.get().get_id() == p.get_id();
            }
        ),
        e.agents.end()
    );

    return;

};

#endif
