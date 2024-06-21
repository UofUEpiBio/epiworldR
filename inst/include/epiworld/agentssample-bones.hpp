#ifndef EPIWORLD_AGENTS_BONES_HPP
#define EPIWORLD_AGENTS_BONES_HPP

class SAMPLETYPE {
public:
    static const int MODEL  = 0;
    static const int ENTITY = 1;
    static const int AGENT  = 2;
};

// template<typename TSeq>
// class Agent;

// template<typename TSeq>
// class Model;

// template<typename TSeq>
// class Entity;

/**
 * @brief Sample of agents
 * 
 * This class allows sampling agents from Entity<TSeq> and Model<TSeq>.
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class AgentsSample {
private:

    size_t sample_size = 0u;

    std::vector< Agent<TSeq>* >* agents = nullptr; ///< Pointer to sample of agents
    size_t * agents_n = nullptr;                    ///< Size of sample of agents
    
    std::vector< size_t >* agents_left = nullptr;  ///< Pointer to agents left (iota)
    size_t * agents_left_n = nullptr;               ///< Size of agents left

    Model<TSeq> * model   = nullptr;   ///< Extracts runif() and (if the case) population.
    Entity<TSeq> * entity = nullptr; ///
    Agent<TSeq> * agent   = nullptr;
    
    int sample_type = SAMPLETYPE::AGENT;
    std::vector< size_t > states = {};

    void sample_n(size_t n); ///< Backbone function for sampling


public:

    // Not available (for now)
    AgentsSample() = delete;                             ///< Default constructor
    AgentsSample(const AgentsSample<TSeq> & a) = delete; ///< Copy constructor
    AgentsSample(AgentsSample<TSeq> && a) = delete;      ///< Move constructor

    AgentsSample(
        Model<TSeq> & model_, size_t n,
        std::vector< size_t > states_ = {},
        bool truncate = false
        );

    AgentsSample(
        Model<TSeq> * model,
        Entity<TSeq> & entity_,
        size_t n,
        std::vector< size_t > states_ = {},
        bool truncate = false
        );

    AgentsSample(
        Model<TSeq> * model,
        Agent<TSeq> & agent_,
        size_t n,
        std::vector< size_t > states_ = {},
        bool truncate = false
        );

    ~AgentsSample();

    typename std::vector< Agent<TSeq> * >::iterator begin();
    typename std::vector< Agent<TSeq> * >::iterator end();

    Agent<TSeq> * operator[](size_t n);
    Agent<TSeq> * operator()(size_t n);
    size_t size() const noexcept;

};

template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> & model_,
    size_t n,
    std::vector< size_t > states_,
    bool truncate
    ) {

    states = states_;

    if (truncate)
    {
        
        if (n > model_.size())
            n = model_.size();

    } else if (n > model_.size())
        throw std::logic_error(
            "There are only " + std::to_string(model_.size()) + " agents. You cannot " +
            "sample " + std::to_string(n));

    sample_size   = n;
    model         = &model_;
    sample_type   = SAMPLETYPE::MODEL;

    agents        = &model_.sampled_population;
    agents_n      = &model_.sampled_population_n;

    agents_left   = &model_.population_left;
    agents_left_n = &model_.population_left_n;

    sample_n(n);
    
    return; 

}

template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> * model,
    Entity<TSeq> & entity_,
    size_t n,
    std::vector< size_t > states_,
    bool truncate
    ) {

    states = states_;

    if (truncate)
    {

        if (n > entity_.size())
            n = entity_.size();

    } else if (n > entity_.size())
        throw std::logic_error(
            "There are only " + std::to_string(entity_.size()) + " agents. You cannot " +
            "sample " + std::to_string(n));

    sample_size   = n;
    model         = &entity_.model;
    sample_type   = SAMPLETYPE::ENTITY;

    agents        = &entity_.sampled_agents;
    agents_n      = &entity_.sampled_agents_n;

    agents_left   = &entity_.sampled_agents_left;
    agents_left_n = &entity_.sampled_agents_left_n;

    sample_n(n);

    return; 

}

/**
 * @brief Sample from the agent's entities
 * 
 * For example, how many individuals the agent contacts in a given point in time.
 * 
 * @tparam TSeq 
 * @param agent_ 
 * @param n Sample size
 * @param truncate If the agent has fewer than `n` connections, then truncate = true
 * will automatically reduce the number of possible samples. Otherwise, if false, then
 * it returns an error.
 */
template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> * model,
    Agent<TSeq> & agent_,
    size_t n,
    std::vector< size_t > states_,
    bool truncate
    )
{

    states = states_;
    sample_type = SAMPLETYPE::AGENT;
    
    agent         = &agent_;

    agents        = &agent_.sampled_agents;
    agents_n      = &agent_.sampled_agents_n;

    // Computing the cumulative sum of counts across entities
    size_t agents_in_entities = 0;
    Entities<TSeq> entities_a = agent->get_entities();

    std::vector< size_t > cum_agents_count(entities_a.size(), 0);
    int idx = -1;
    for (auto & e : entities_a)
    {
        if (++idx == 0)
            cum_agents_count[idx] = (e->size() - 1u);
        else
            cum_agents_count[idx] = (
                (e->size() - 1u) + 
                cum_agents_count[idx - 1]
            );

        agents_in_entities += (e->size() - 1u);
    }

    if (truncate)
    {
        
        if (n > agents_in_entities)
            n = agents_in_entities;

    } else if (n > agents_in_entities)
        throw std::logic_error(
            "There are only " + std::to_string(agents_in_entities) +
            " agents. You cannot " +
            "sample " + std::to_string(n)
            );

    sample_size = n;

    if (agents->size() < n)
        agents->resize(n);

    size_t i_obs = 0u;
    for (size_t i = 0u; i < sample_size; ++i)
    {

        // Sampling a single agent from the set of entities
        int jth = std::floor(model->runif() * agents_in_entities);
        for (size_t e = 0u; e < cum_agents_count.size(); ++e)
        {
            
            // Are we in the limit?
            if (jth <= cum_agents_count[e])
            {
                size_t agent_idx = 0u;
                if (e == 0) // From the first group
                    agent_idx = entities_a[e][jth]->get_id();
                else
                    agent_idx = entities_a[e][jth - cum_agents_count[e - 1]]->get_id();


                // Checking if states was specified
                if (states.size())
                {

                    // Getting the state
                    size_t state = model->population[agent_idx].get_state();

                    if (std::find(states.begin(), states.end(), state) != states.end())
                        continue;

                }
                
                agents->operator[](i_obs++) = &(model->population[agent_idx]);

                break;
            }

        }
    }

    return; 

}

template<typename TSeq>
inline AgentsSample<TSeq>::~AgentsSample() {}

template<typename TSeq>
inline size_t AgentsSample<TSeq>::size() const noexcept
{ 
    return this->sample_size;
}

template<typename TSeq>
inline Agent<TSeq> * AgentsSample<TSeq>::operator[](size_t i)
{

    return agents->operator[](i);

}

template<typename TSeq>
inline Agent<TSeq> * AgentsSample<TSeq>::operator()(size_t i)
{

    if (i >= this->sample_size)
        throw std::range_error("The requested agent is out of range.");

    return agents->operator[](i);

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator AgentsSample<TSeq>::begin()
{

    if (sample_size > 0u)
        return agents->begin();
    else
        return agents->end();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator AgentsSample<TSeq>::end()
{

    return agents->begin() + sample_size;

}

template<typename TSeq>
inline void AgentsSample<TSeq>::sample_n(size_t n)
{

    // Reducing size
    if (states.size())
    {
            
        // Getting the number of agents left
        agents_left->clear();

        if (sample_type == SAMPLETYPE::MODEL)
        {

            // Making some room
            agents_left->reserve(model->size());

            // Iterating through the agents in the population
            for (size_t a_i = 0u; a_i < model->population.size(); ++a_i)
            {

                // If the agent is within the selected set of states,
                // then we add it to the list of agents left
                size_t s = model->population[a_i].get_state();
                if (std::find(states.begin(), states.end(), s) != states.end())
                    agents_left->push_back(a_i);

            }

        } 
        else if (sample_type == SAMPLETYPE::ENTITY)
        {

            // Making room
            agents_left->reserve(entity->size());

            // Iterating through the agents in the entity
            for (size_t a_i = 0u; a_i < entity->size(); ++a_i)
            {
                size_t s = model->population[entity->agents[a_i]].get_state();
                if (std::find(states.begin(), states.end(), s) != states.end())
                    agents_left->push_back(a_i);

            }

        }

    } else {

        // Checking if the size of the entity has changed (or hasn't been initialized)
        if (sample_type == SAMPLETYPE::MODEL)
        {

            if (model->size() != agents_left->size())
            {
                agents_left->resize(model->size(), 0u);
                std::iota(agents_left->begin(), agents_left->end(), 0u);
            }

        } else if (sample_type == SAMPLETYPE::ENTITY) {

            if (entity->size() != agents_left->size())
            {

                agents_left->resize(entity->size(), 0u);
                std::iota(agents_left->begin(), agents_left->end(), 0u);

            }

        } 

    }

    // Restart the counter of agents left
    *agents_left_n = agents_left->size();

    // Making sure we have enough room for the sample of agents
    if (agents->size() < sample_size)
        agents->resize(sample_size, nullptr);

    if (sample_type == SAMPLETYPE::MODEL)
    {

        #ifdef EPI_DEBUG
        std::vector< bool > __sampled(model->size(), true);
        for (auto & a_i: *agents_left)
            __sampled[a_i] = false;
        #endif

        for (size_t i = 0u; i < n; ++i)
        {

            // Sampling from 0 to (agents_left_n - 1)
            size_t ith_ = static_cast<size_t>(model->runif() * ((*agents_left_n)--));
            
            // Getting the id of the agent and adding it to the list of agents
            size_t ith  = agents_left->operator[](ith_);
            agents->operator[](i) = &model->population[ith];

            #ifdef EPI_DEBUG
            if (__sampled[ith])
                throw std::logic_error("The same agent was sampled twice.");
            else
                __sampled[ith] = true;
            #endif

            // Updating list
            std::swap(
                agents_left->operator[](ith_),
                agents_left->operator[](*agents_left_n)
                );

        }


    }
    else if (sample_type == SAMPLETYPE::ENTITY)
    {

        #ifdef EPI_DEBUG
        std::vector< bool > __sampled(entity->size(), true);
        for (auto & a_i: *agents_left)
            __sampled[a_i] = false;
        #endif

        for (size_t i = 0u; i < n; ++i)
        {

            size_t ith_ = static_cast<size_t>(model->runif() * ((*agents_left_n)--));
            size_t ith  = agents_left->operator[](ith_);
            agents->operator[](i) = &model->population[entity->agents[ith]];

            #ifdef EPI_DEBUG
            if (__sampled[ith])
                throw std::logic_error("The same agent was sampled twice.");
            else
                __sampled[ith] = true;
            #endif

            // Updating list
            std::swap(agents_left->operator[](ith_), agents_left->operator[](*agents_left_n));

        }

    }

    return;

}

#endif