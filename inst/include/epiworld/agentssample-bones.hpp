#ifndef EPIWORLD_AGENTS_BONES_HPP
#define EPIWORLD_AGENTS_BONES_HPP

/**
 * @brief Sample of agents
 * 
 * A sample without replacement from the model's agents, optionally restricted
 * to some states.
 * 
 * @tparam TSeq 
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class AgentsSample {
private:

    size_t sample_size = 0u;

    std::vector< Agent<TSeq>* >* agents = nullptr; ///< Pointer to sample of agents
    size_t * agents_n = nullptr;                    ///< Size of sample of agents
    
    std::vector< size_t >* agents_left = nullptr;  ///< Pointer to agents left (iota)
    size_t * agents_left_n = nullptr;               ///< Size of agents left

    Model<TSeq> * model = nullptr; ///< Extracts runif() and the population.
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

    agents        = &model_.sampled_population;
    agents_n      = &model_.sampled_population_n;

    agents_left   = &model_.population_left;
    agents_left_n = &model_.population_left_n;

    sample_n(n);
    
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
    // Check for null pointer
    if (agents == nullptr)
        return typename std::vector< Agent<TSeq> * >::iterator{};
    
    return agents->begin();
}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator AgentsSample<TSeq>::end()
{
    // Check for null pointer
    if (agents == nullptr)
        return typename std::vector< Agent<TSeq> * >::iterator{};
    
    // Ensure we don't go beyond the actual vector size
    size_t actual_end = std::min(sample_size, agents->size());
    return agents->begin() + actual_end;
}

template<typename TSeq>
inline void AgentsSample<TSeq>::sample_n(size_t n)
{

    if (states.size())
    {

        // Only the agents in the selected states can be sampled
        agents_left->clear();
        agents_left->reserve(model->size());
        for (size_t a_i = 0u; a_i < model->population.size(); ++a_i)
        {
            size_t s = model->population[a_i].get_state();
            if (std::find(states.begin(), states.end(), s) != states.end())
                agents_left->push_back(a_i);
        }

    }
    else if (model->size() != agents_left->size())
    {

        // The population changed size (or this is the first sample)
        agents_left->resize(model->size(), 0u);
        std::iota(agents_left->begin(), agents_left->end(), 0u);

    }

    // Restart the counter of agents left
    *agents_left_n = agents_left->size();

    // Making sure we have enough room for the sample of agents
    if (agents->size() < sample_size)
        agents->resize(sample_size, nullptr);

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

    return;

}

#endif
