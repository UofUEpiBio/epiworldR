#ifndef EPIWORLD_AGENT_EVENTS_MEAT_HPP
#define EPIWORLD_AGENT_EVENTS_MEAT_HPP

template<typename TSeq>
inline void Model<TSeq>::_event_add_virus(Event<TSeq> & a)
{

    Agent<TSeq> *  p = a.agent;
    VirusPtr<TSeq> & v = a.virus;
    
    db.record_transmission(
        v->get_agent() ? v->get_agent()->get_id() : -1,
        p->get_id(),
        v->get_id(),
        v->get_date() 
    );
    
    p->virus = std::move(v);
    p->virus->set_date(today());
    p->virus->set_agent(p);

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
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
    db.today_virus.at(p->virus->get_id()).at(
        a.new_state != -99 ? a.new_state : p->state
    )++;
    #else
    db.today_virus[p->virus->get_id()][
        a.new_state != -99 ? a.new_state : p->state
    ]++;
    #endif

}

template<typename TSeq>
inline void Model<TSeq>::_event_add_tool(Event<TSeq> & a)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> & t = a.tool;
    
    // Update tool accounting
    p->tools.emplace_back(std::move(t));
    size_t tool_pos = p->tools.size() - 1u;

    p->tools[tool_pos]->set_date(today());
    p->tools[tool_pos]->set_agent(p, tool_pos);

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && static_cast<int>(p->state) != a.new_state)
    {
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

    db.today_tool[p->tools[tool_pos]->get_id()][
        a.new_state != -99 ? a.new_state : p->state
    ]++;


}

template<typename TSeq>
inline void Model<TSeq>::_event_rm_virus(Event<TSeq> & a)
{

    Agent<TSeq> * p    = a.agent;
    VirusPtr<TSeq> & v = a.virus;

    // Calling the virus action over the removed virus
    v->post_recovery(this);

    p->virus = nullptr;

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
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

    // The counter of the removed virus must be decremented from the
    // agent's current state (pre-transition), not state_prev.
    // state_prev may refer to an earlier same-day state.
    #ifdef EPI_DEBUG
    db.today_virus.at(v->get_id()).at(p->state)--;
    #else
    db.today_virus[v->get_id()][p->state]--;
    #endif

    
    return;

}

template<typename TSeq>
inline void Model<TSeq>::_event_rm_tool(Event<TSeq> & a)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> & t = a.tool;
    bool removed = false;

    if (t)
    {
        // Remove tool(s) by id, following an erase/remove approach.
        auto new_end = std::remove_if(
            p->tools.begin(),
            p->tools.end(),
            [&t](const ToolPtr<TSeq> & tool_ptr) -> bool
            {
                return tool_ptr && (tool_ptr->get_id() == t->get_id());
            }
        );

        removed = (new_end != p->tools.end());
        p->tools.erase(new_end, p->tools.end());

        for (size_t i = 0u; i < p->tools.size(); ++i)
            p->tools[i]->pos_in_agent = static_cast<int>(i);
    }

    // Change of state needs to be recorded and updated on the
    // tools.
    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
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

    // Like rm_virus, remove the tool count from the agent's current
    // state (pre-transition). Using state_prev can underflow when
    // the agent changed state earlier in the day.
    if (removed)
    {
        #ifdef EPI_DEBUG
        db.today_tool.at(t->get_id()).at(p->state)--;
        #else
        db.today_tool[t->get_id()][p->state]--;
        #endif

        t->agent = nullptr;
        t->pos_in_agent = -99;
    }

    return;

}

template<typename TSeq>
inline void Model<TSeq>::_event_change_state(Event<TSeq> & a)
{

    Agent<TSeq> * p = a.agent;

    if ((a.new_state != -99) && (static_cast<int>(p->state) != a.new_state))
    {
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
inline void Model<TSeq>::_event_add_entity(Event<TSeq> & a)
{

    Agent<TSeq> *  p = a.agent;
    Entity<TSeq> * e = a.entity;

    // Checking the agent and the entity are not linked
    if ((p->get_n_entities() > 0) && (e->size() > 0))
    {

        if (p->get_n_entities() > e->size()) // Slower search through the agent
        {
            for (size_t agent_id : e->get_agents())
                if(static_cast<int>(agent_id) == p->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }
        else                                 // Slower search through the entity
        {
            for (size_t entity_id : p->get_entities())
                if(static_cast<int>(entity_id) == e->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }

        // It means that agent and entity were not associated.
    }

    // Adding the to agent and the entity
    p->entities.push_back(static_cast<size_t>(e->get_id()));
    e->agents.push_back(static_cast<size_t>(p->get_id()));

}

template<typename TSeq>
inline void Model<TSeq>::_event_rm_entity(Event<TSeq> & a)
{

    Agent<TSeq> &  p = *a.agent;
    Entity<TSeq> & e = *a.entity;

    // Remove entity from agent's entity list
    p.entities.erase(
        std::remove(
            p.entities.begin(),
            p.entities.end(),
            static_cast<size_t>(e.get_id())
        ),
        p.entities.end()
    );

    // Remove agent from entity's agent list
    e.agents.erase(
        std::remove(
            e.agents.begin(),
            e.agents.end(),
            static_cast<size_t>(p.get_id())
        ),
        e.agents.end()
    );

    return;

}

#endif
