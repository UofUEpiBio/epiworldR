#ifndef EPIWORLD_ENTITY_MEAT_HPP
#define EPIWORLD_ENTITY_MEAT_HPP

template <typename TSeq = EPI_DEFAULT_TSEQ>
inline EntityToAgentFun<TSeq> entity_to_unassigned_agents()
{

    return [](Entity<TSeq> & e, Model<TSeq> * m) -> void {

        
        // Preparing the sampling space
        std::vector< size_t > idx;
        for (const auto & a: m->get_agents())
            if (a.get_n_entities() == 0)
                idx.push_back(a.get_id());
        size_t n = idx.size();

        // Figuring out how many to sample
        int n_to_sample;
        if (e.get_prevalence_as_proportion())
        {
            n_to_sample = static_cast<int>(std::floor(e.get_prevalence() * n));
            if (n_to_sample > static_cast<int>(n))
                --n_to_sample;

        } else
        {
            n_to_sample = static_cast<int>(e.get_prevalence());
            if (n_to_sample > static_cast<int>(n))
                throw std::range_error("There are only " + std::to_string(n) + 
                " individuals in the population. Cannot add the entity to " +
                    std::to_string(n_to_sample));
        }

        int n_left = n;
        for (int i = 0; i < n_to_sample; ++i)
        {
            int loc = static_cast<epiworld_fast_uint>(
                floor(m->runif() * n_left--)
                );

            // Correcting for possible overflow
            if ((loc > 0) && (loc >= n_left))
                loc = n_left - 1;

            m->get_agent(idx[loc]).add_entity(e, m);

            std::swap(idx[loc], idx[n_left]);

        }

    };

}

template<typename TSeq = int>
inline EntityToAgentFun<TSeq> entity_to_agent_range(
    int from,
    int to,
    bool to_unassigned = false
    ) {

    if (to_unassigned)
    {

        return [from, to](Entity<TSeq> & e, Model<TSeq> * m) -> void {

            auto & agents = m->get_agents();
            for (size_t i = from; i < to; ++i)
            {
                if (agents[i].get_n_entities() == 0)
                    e.add_agent(&agents[i], m);
                else
                    throw std::logic_error(
                        "Agent " + std::to_string(i) + " already has an entity."
                    );
            }
            
            return;

        };

    }
    else
    {

        return [from, to](Entity<TSeq> & e, Model<TSeq> * m) -> void {

            auto & agents = m->get_agents();
            for (size_t i = from; i < to; ++i)
            {
                e.add_agent(&agents[i], m);
            }
            
            return;

        };

    }
}

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
inline void Entity<TSeq>::rm_agent(size_t idx)
{
    if (idx >= n_agents)
        throw std::out_of_range(
            "Trying to remove agent "+ std::to_string(idx) +
            " out of " + std::to_string(n_agents)
            );

    model->get_agents()[agents[idx]].rm_entity(*this, model);

    return;
}

template<typename TSeq>
inline size_t Entity<TSeq>::size() const noexcept
{
    return n_agents;
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
inline typename std::vector< Agent<TSeq> * >::iterator Entity<TSeq>::begin()
{

    if (n_agents == 0)
        return agents.end();

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator Entity<TSeq>::end()
{
    return agents.begin() + n_agents;
}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::const_iterator Entity<TSeq>::begin() const
{

    if (n_agents == 0)
        return agents.end();

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::const_iterator Entity<TSeq>::end() const
{
    return agents.begin() + n_agents;
}

template<typename TSeq>
inline Agent<TSeq> * Entity<TSeq>::operator[](size_t i)
{
    if (n_agents <= i)
        throw std::logic_error(
            "There are not that many agents in this entity. " +
            std::to_string(n_agents) + " <= " + std::to_string(i)
            );

    return &model->get_agents()[i];
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
    sampled_agents.clear();
    sampled_agents_n = 0u;
    sampled_agents_left.clear();
    sampled_agents_left_n = 0u;

    this->agents.clear();
    this->n_agents = 0u;
    this->agents_location.clear();

    return;

}

template<typename TSeq>
inline bool Entity<TSeq>::operator==(const Entity<TSeq> & other) const
{

    if (id != other.id)
        return false;

    if (n_agents != other.n_agents)
        return false;

    for (size_t i = 0u; i < n_agents; ++i)
    {
        if (agents[i] != other.agents[i])
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
inline void Entity<TSeq>::distribute()
{

    // Starting first infection
    int n = this->model->size();
    std::vector< size_t > idx(n);

    if (dist_fun)
    {

        dist_fun(*this, model);

    }
    else
    {

        // Picking how many
        int n_to_assign;
        if (prevalence_as_proportion)
        {
            n_to_assign = static_cast<int>(std::floor(prevalence * size()));
        }
        else
        {
            n_to_assign = static_cast<int>(prevalence);
        }

        if (n_to_assign > static_cast<int>(model->size()))
            throw std::range_error("There are only " + std::to_string(model->size()) + 
            " individuals in the population. Cannot add the entity to " + std::to_string(n_to_assign));
        
        int n_left = n;
        std::iota(idx.begin(), idx.end(), 0);
        while ((n_to_assign > 0) && (n_left > 0))
        {
            int loc = static_cast<epiworld_fast_uint>(
                floor(model->runif() * n_left--)
                );

            // Correcting for possible overflow
            if ((loc > 0) && (loc >= n_left))
                loc = n_left - 1;
            
            auto & agent = model->get_agent(idx[loc]);

            if (!agent.has_entity(id))
            {
                agent.add_entity(
                    *this, this->model, this->state_init, this->queue_init
                    );
                n_to_assign--;
            }
                            
            std::swap(idx[loc], idx[n_left]);

        }

    }

}

template<typename TSeq>
inline std::vector< size_t > & Entity<TSeq>::get_agents()
{
    return agents;
}

template<typename TSeq>
inline void Entity<TSeq>::print() const
{

    printf_epiworld(
        "Entity '%s' (id %i) with %i agents.\n",
        this->entity_name.c_str(),
        static_cast<int>(id),
        static_cast<int>(n_agents)
    );
}

template<typename TSeq>
inline void Entity<TSeq>::set_prevalence(
    epiworld_double p,
    bool as_proportion
)
{
    prevalence = p;
    prevalence_as_proportion = as_proportion;
}

template<typename TSeq>
inline epiworld_double Entity<TSeq>::get_prevalence() const noexcept
{
    return prevalence;
}

template<typename TSeq>
inline bool Entity<TSeq>::get_prevalence_as_proportion() const noexcept
{
    return prevalence_as_proportion;
}

template<typename TSeq>
inline void Entity<TSeq>::set_dist_fun(EntityToAgentFun<TSeq> fun)
{
    dist_fun = fun;
}

#endif