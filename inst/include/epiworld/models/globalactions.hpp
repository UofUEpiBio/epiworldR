#ifndef EPIWORLD_GLOBALACTIONS_HPP
#define EPIWORLD_GLOBALACTIONS_HPP

// This function creates a global action that distributes a tool
// to agents with probability p.
/**
 * @brief Global action that distributes a tool to agents with probability p.
 * 
 * @tparam TSeq Sequence type (should match `TSeq` across the model)
 * @param p Probability of distributing the tool.
 * @param tool Tool function.
 * @return std::function<void(Model<TSeq>*)> 
 */
template<typename TSeq>
inline std::function<void(Model<TSeq>*)> globalaction_tool(
    Tool<TSeq> & tool,
    double p
) {

    std::function<void(Model<TSeq>*)> fun = [p,&tool](
        Model<TSeq> * model
        ) -> void {

        for (auto & agent : model->get_agents())
        {

            // Check if the agent has the tool
            if (agent.has_tool(tool))
                continue;

            // Adding the tool
            if (model->runif() < p)
                agent.add_tool(tool, model);
            
        
        }

        #ifdef EPIWORLD_DEBUG
        tool.print();
        #endif

        return;
            

    };

    return fun;

}

// Same function as above, but p is now a function of a vector of coefficients
// and a vector of variables.
/**
 * @brief Global action that distributes a tool to agents with probability
 * p = 1 / (1 + exp(-\sum_i coef_i * agent(vars_i))).
 * 
 * @tparam TSeq Sequence type (should match `TSeq` across the model)
 * @param coefs Vector of coefficients.
 * @param vars Vector of variables.
 * @param tool_fun Tool function.
 * @return std::function<void(Model<TSeq>*)> 
 */
template<typename TSeq>
inline std::function<void(Model<TSeq>*)> globalaction_tool_logit(
    Tool<TSeq> & tool,
    std::vector< size_t > vars,
    std::vector< double > coefs
) {

    std::function<void(Model<TSeq>*)> fun = [coefs,vars,&tool](
        Model<TSeq> * model
        ) -> void {

        for (auto & agent : model->get_agents())
        {

            // Check if the agent has the tool
            if (agent.has_tool(tool))
                continue;

            // Computing the probability using a logit. Uses OpenMP reduction
            // to sum the coefficients.
            double p = 0.0;
            #pragma omp parallel for reduction(+:p)
            for (size_t i = 0u; i < coefs.size(); ++i)
                p += coefs.at(i) * agent(vars[i]);

            p = 1.0 / (1.0 + std::exp(-p));

            // Adding the tool
            if (model->runif() < p)
                agent.add_tool(tool, model);
            
        
        }

        #ifdef EPIWORLD_DEBUG
        tool.print();
        #endif

        return;
            

    };

    return fun;

}

// A global action that updates a parameter in the model.
/**
 * @brief Global action that updates a parameter in the model.
 * 
 * @tparam TSeq Sequence type (should match `TSeq` across the model)
 * @param param Parameter to update.
 * @param value Value to update the parameter to.
 * @return std::function<void(Model<TSeq>*)> 
 */
template<typename TSeq>
inline std::function<void(Model<TSeq>*)> globalaction_set_param(
    std::string param,
    double value
) {

    std::function<void(Model<TSeq>*)> fun = [value,param](
        Model<TSeq> * model
        ) -> void {

        model->set_param(param, value);

        return;
            

    };

    return fun;

}
#endif