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
            if ((loc > 0) && (loc >= n_available))
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

#endif