#ifndef EPIWORLD_VIRUS_DISTRIBUTE_MEAT_HPP
#define EPIWORLD_VIRUS_DISTRIBUTE_MEAT_HPP

/**
 * Distributes a virus to a set of agents.
 *
 * This function takes a vector of agent IDs and returns a lambda function that
 * can be used to distribute a virus to the specified agents.
 *
 * @param agents_ids A vector of agent IDs representing the set of agents to
 * distribute the virus to.
 * 
 * @return A lambda function that takes a Virus object and a Model object and
 * distributes the virus to the specified agents.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline VirusToAgentFun<TSeq> distribute_virus_to_set(
    std::vector< size_t > agents_ids
) {

    return [agents_ids](
        Virus<TSeq> & virus, Model<TSeq> * model
    ) -> void 
    { 
        // Adding action
        for (auto i: agents_ids)
        {
            model->get_agent(i).set_virus(
                virus,
                const_cast<Model<TSeq> * >(model)
                );
        }
    };

}

/**
 * @brief Distributes a virus randomly to agents.
 * 
 * This function takes a sequence of agents and randomly assigns a virus to
 * each agent.
 * 
 * @tparam TSeq The type of the sequence of agents.
 * @param agents The sequence of agents to distribute the virus to.
 * @return A function object that assigns a virus to each agent randomly.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline VirusToAgentFun<TSeq> distribute_virus_randomly(
    epiworld_double prevalence,
    bool prevalence_as_proportion = true,
    std::vector< size_t > agents_ids = {}
) {

    auto agents_ids_ptr = std::make_shared< std::vector< size_t > >(agents_ids);

    return [prevalence,prevalence_as_proportion,agents_ids_ptr](
        Virus<TSeq> & virus, Model<TSeq> * model
    ) -> void 
    { 
        
        // Figuring out how what agents are available
        bool use_set = agents_ids_ptr->size() > 0;
        std::vector< size_t > idx;
        if (use_set)
        {
            for (const auto & agent: *agents_ids_ptr)
                if (model->get_agent(agent).get_virus() == nullptr)
                    idx.push_back(agent);
        }
        else 
        {
            for (const auto & agent: model->get_agents())
                if (agent.get_virus() == nullptr)
                    idx.push_back(agent.get_id());
        }

        // Picking how many
        int n = use_set ?
            static_cast<int>(idx.size()) :
            static_cast<int>(model->size())
            ;
            
        int n_available = static_cast<int>(idx.size());
        int n_to_sample;
        if (prevalence_as_proportion)
        {
            n_to_sample = static_cast<int>(std::floor(
                prevalence * static_cast< epiworld_double >(n)
            ));

            // Correcting for possible overflow
            if (n_to_sample > n)
                n_to_sample = n;
        }
        else
        {
            n_to_sample = static_cast<int>(prevalence);
        }

        if (n_to_sample > n_available)
            throw std::range_error(
                "There are only " + std::to_string(n_available) + 
                " individuals with no virus in the population. " +
                "Cannot add the virus to " +
                std::to_string(n_to_sample)
            );
        
        auto & population = model->get_agents();
        for (int i = 0; i < n_to_sample; ++i)
        {

            int loc = static_cast<epiworld_fast_uint>(
                floor(model->runif() * (n_available--))
                );

            // Correcting for possible overflow
            if ((n_available > 0) && (loc >= n_available))
                loc = n_available - 1;

            Agent<TSeq> & agent = population[idx[loc]];
            
            // Adding action
            agent.set_virus(
                virus,
                const_cast<Model<TSeq> * >(model)
                );

            // Adjusting sample
            std::swap(idx[loc], idx[n_available]);

        }

    };

}

/**
 * Function template to distribute a virus to agents in each entity of a model.
 * 
 * @tparam TSeq The sequence type used in the model.
 * @param prevalence A vector of prevalences for each entity in the model.
 * @param as_proportion Flag indicating whether the prevalences are given as
 * proportions or absolute values.
 * @return A lambda function that distributes the virus to agents in each entity
 * of the model.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline VirusToAgentFun<TSeq> distribute_virus_to_entities(
    std::vector< double > prevalence,
    bool as_proportion = true
)
{

    // Checking proportions are within range
    if (prevalence.size() == 0)
        throw std::range_error("Prevalence vector cannot be empty");

    if (as_proportion)
    {
        for (auto p: prevalence)
        {
            if ((p < 0.0) || (p > 1.0))
                throw std::range_error("Proportions must be between 0 and 1");
        }
    }
    else
    {
        for (auto p: prevalence)
        {
            if (p < 0.0)
                throw std::range_error("Count values in prevalence must be non-negative");
        }
    }

    return [prevalence, as_proportion](
        Virus<TSeq> & virus, Model<TSeq> * model
    ) -> void 
    { 

        // Checking the number of entities
        if (prevalence.size() != model->get_entities().size())
            throw std::range_error(
                "Prevalence vector size (" + std::to_string(prevalence.size()) +
                ") does not match number of entities (" + std::to_string(model->get_entities().size()) + ")"
                );

        // Adding action
        auto & entities = model->get_entities();
        auto & population = model->get_agents();
        for (size_t e = 0; e < entities.size(); ++e)
        {
            auto & entity_e = entities[e];
            auto & agent_ids = entity_e.get_agents();
            double prevalence_e = prevalence[e];
            size_t n = agent_ids.size();
            size_t n_to_distribute;
            if (as_proportion)
                n_to_distribute = static_cast<size_t>(std::floor(prevalence_e * n));
            else
                n_to_distribute = static_cast<size_t>(prevalence_e);

            if (n_to_distribute > n)
                n_to_distribute = n;

            std::vector< size_t > idx = agent_ids;
            for (size_t i = 0u; i < n_to_distribute; ++i)
            {
                size_t loc = static_cast<epiworld_fast_uint>(
                    floor(model->runif() * n--)
                    );

                if ((n > 0) && (loc >= n))
                    loc = n - 1;
                
                population[idx[loc]].set_virus(
                    virus,
                    const_cast< Model<TSeq> * >(model)
                    );
                
                std::swap(idx[loc], idx[n]);

            }
        }
    };

}


#endif