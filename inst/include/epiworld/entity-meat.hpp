#ifndef EPIWORLD_ENTITY_MEAT_HPP
#define EPIWORLD_ENTITY_MEAT_HPP

template<typename TSeq>
inline void Entity<TSeq>::add_agent(
    Agent<TSeq> & p,
    Model<TSeq> * model
    )
{

    // Need to add it to the events, through the individual
    p.add_entity(*this, model);    

}

template<typename TSeq>
inline void Entity<TSeq>::add_agent(
    Agent<TSeq> * p,
    Model<TSeq> * model
    )
{
    p->add_entity(*this, model);
}

template<typename TSeq>
inline void Entity<TSeq>::rm_agent(size_t idx, Model<TSeq> * model)
{
    if (idx >= size())
        throw std::out_of_range(
            "Trying to remove agent "+ std::to_string(idx) +
            " out of " + std::to_string(size())
            );

    agents[idx].rm_entity(*this, model);

    return;
}

template<typename TSeq>
inline size_t Entity<TSeq>::size() const noexcept
{
    return agents.size();
}

template<typename TSeq>
inline void Entity<TSeq>::set_location(std::vector< epiworld_double > loc)
{
    location = loc;
}

template<typename TSeq>
inline std::vector< epiworld_double > & Entity<TSeq>::get_location()
{
    return location;
}

template<typename TSeq>
inline typename std::vector< std::reference_wrapper<Agent<TSeq>> >::iterator Entity<TSeq>::begin()
{

    if (agents.size() == 0)
        return typename std::vector< std::reference_wrapper<Agent<TSeq>> >::iterator{};

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< std::reference_wrapper<Agent<TSeq>> >::iterator Entity<TSeq>::end()
{
    if (agents.size() == 0)
        return typename std::vector< std::reference_wrapper<Agent<TSeq>> >::iterator{};

    return agents.begin() + agents.size();
}

template<typename TSeq>
inline typename std::vector< std::reference_wrapper<Agent<TSeq>> >::const_iterator Entity<TSeq>::begin() const
{

    if (agents.size() == 0)
        return typename std::vector< std::reference_wrapper<Agent<TSeq>> >::const_iterator{};

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< std::reference_wrapper<Agent<TSeq>> >::const_iterator Entity<TSeq>::end() const
{
    if (agents.size() == 0)
        return typename std::vector< std::reference_wrapper<Agent<TSeq>> >::const_iterator{};

    return agents.begin() + agents.size();
}

template<typename TSeq>
size_t Entity<TSeq>::operator[](size_t i)
{
    if (agents.size() <= i)
        throw std::logic_error(
            "There are not that many agents in this entity. " +
            std::to_string(agents.size()) + " <= " + std::to_string(i)
            );

    return i;
}

template<typename TSeq>
inline int Entity<TSeq>::get_id() const noexcept
{
    return id;
}

template<typename TSeq>
inline const std::string & Entity<TSeq>::get_name() const noexcept
{
    return entity_name;
}

template<typename TSeq>
inline void Entity<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    state_init = init;
    state_post = end;
}

template<typename TSeq>
inline void Entity<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    queue_init = init;
    queue_post = end;
}

template<typename TSeq>
inline void Entity<TSeq>::get_state(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = state_init;

    if (post != nullptr)
        *post = state_post;

}

template<typename TSeq>
inline void Entity<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = queue_init;

    if (post != nullptr)
        *post = queue_post;

}

template<typename TSeq>
inline void Entity<TSeq>::reset()
{

    this->agents.clear();

    return;

}

template<typename TSeq>
inline bool Entity<TSeq>::operator==(const Entity<TSeq> & other) const
{

    if (id != other.id)
        return false;

    if (agents.size() != other.agents.size())
        return false;

    for (size_t i = 0u; i < agents.size(); ++i)
    {
        if (agents[i].get() != other.agents[i].get())
            return false;
    }


    if (max_capacity != other.max_capacity)
        return false;

    if (entity_name != other.entity_name)
        return false;

    if (location.size() != other.location.size())
        return false;

    for (size_t i = 0u; i < location.size(); ++i)
    {

        if (location[i] != other.location[i])
            return false;

    }

    if (state_init != other.state_init)
        return false;

    if (state_post != other.state_post)
        return false;

    if (queue_init != other.queue_init)
        return false;

    if (queue_post != other.queue_post)
        return false;

    return true;

}

template<typename TSeq>
inline void Entity<TSeq>::distribute(Model<TSeq> * model)
{

    if (dist_fun)
    {

        dist_fun(*this, model);

    }

}

template<typename TSeq>
inline std::vector< std::reference_wrapper<Agent<TSeq>> > & Entity<TSeq>::get_agents()
{
    return agents;
}

template<typename TSeq>
inline std::vector< size_t > Entity<TSeq>::get_agents_ids() const
{
    std::vector< size_t > res;
    for (const Agent<TSeq> & agent: agents)
        res.push_back(agent.get_id());

    return res;
}

template<typename TSeq>
inline void Entity<TSeq>::print() const
{

    printf_epiworld(
        "Entity '%s' (id %i) with %i agents.\n",
        this->entity_name.c_str(),
        static_cast<int>(id),
        static_cast<int>(agents.size())
    );
}

template<typename TSeq>
inline void Entity<TSeq>::set_distribution(EntityToAgentFun<TSeq> fun)
{
    dist_fun = fun;
}

#endif