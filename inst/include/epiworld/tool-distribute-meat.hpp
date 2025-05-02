#ifndef TOOL_DISTRIBUTE_MEAT_HPP
#define TOOL_DISTRIBUTE_MEAT_HPP

/**
 * @brief Distributes a tool to a set of agents.
 * 
 * This function takes a vector of agent IDs and returns a lambda function that
 * distributes a tool to each agent in the set.
 * 
 * The lambda function takes a reference to a Tool object and a pointer to a
 * Model object as parameters. It iterates over the agent IDs and adds the tool
 * to each agent using the add_tool() method of the Model object.
 * 
 * @tparam TSeq The sequence type used in the Tool and Model objects.
 * @param agents_ids A vector of agent IDs representing the set of agents to
 * distribute the tool to.
 * @return A lambda function that distributes the tool to the set of agents.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline ToolToAgentFun<TSeq> distribute_tool_to_set(
    std::vector< size_t > agents_ids
) {

    return [agents_ids](
        Tool<TSeq> & tool, Model<TSeq> * model
    ) -> void 
    { 
        // Adding action
        for (auto i: agents_ids)
        {
            model->get_agent(i).add_tool(
                tool,
                const_cast<Model<TSeq> * >(model)
                );
        }
    };

}

/**
 * Function template to distribute a tool randomly to agents in a model.
 * 
 * @tparam TSeq The sequence type used in the model.
 * @param prevalence The prevalence of the tool in the population.
 * @param as_proportion Flag indicating whether the prevalence is given as a
 * proportion or an absolute value.
 * @return A lambda function that distributes the tool randomly to agents in
 * the model.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline ToolToAgentFun<TSeq> distribute_tool_randomly(
    epiworld_double prevalence,
    bool as_proportion = true,
    std::vector< size_t > agents_ids = {}
) {

    auto agents_ids_ptr = std::make_shared< std::vector< size_t > >(agents_ids);

    return [prevalence,as_proportion,agents_ids_ptr](
        Tool<TSeq> & tool, Model<TSeq> * model
        ) -> void {

            // Figuring out how what agents are available
            bool use_set = agents_ids_ptr->size() > 0;

            // Picking how many
            int n_to_distribute;
            int n = use_set ? 
                static_cast<int>(agents_ids_ptr->size()) :
                static_cast<int>(model->size());
                
            if (as_proportion)
            {
                n_to_distribute = static_cast<int>(std::floor(prevalence * n));

                // Correcting for possible rounding errors
                if (n_to_distribute > n)
                    n_to_distribute = n;

            }
            else
            {
                n_to_distribute = static_cast<int>(prevalence);
            }

            if (n_to_distribute > n)
                throw std::range_error("There are only " + std::to_string(n) + 
                " individuals in the population. Cannot add the tool to " + std::to_string(n_to_distribute));
            
            std::vector< int > idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            auto & population = model->get_agents();
            for (int i = 0u; i < n_to_distribute; ++i)
            {
                int loc = static_cast<epiworld_fast_uint>(
                    floor(model->runif() * n--)
                    );

                if ((loc > 0) && (loc == n))
                    loc--;
                
                population[idx[loc]].add_tool(
                    tool,
                    const_cast< Model<TSeq> * >(model)
                    );
                
                std::swap(idx[loc], idx[n]);

            }

        };

}
#endif