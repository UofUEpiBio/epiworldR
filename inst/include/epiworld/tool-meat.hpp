
#ifndef EPIWORLD_TOOLS_MEAT_HPP
#define EPIWORLD_TOOLS_MEAT_HPP

/**
 * @brief Factory function of ToolFun base on logit
 * 
 * @tparam TSeq 
 * @param vars Vector indicating the position of the variables to use.
 * @param coefs Vector of coefficients.
 * @return ToolFun<TSeq> 
 */
template<typename TSeq>
inline ToolFun<TSeq> tool_fun_logit(
    std::vector< int > vars,
    std::vector< double > coefs,
    Model<TSeq> * model
) {

    // Checking that there are features
    if (coefs.size() == 0u)
        throw std::logic_error(
            "The -coefs- argument should feature at least one element."
            );

    if (coefs.size() != vars.size())
        throw std::length_error(
            std::string("The length of -coef- (") +
            std::to_string(coefs.size()) + 
            std::string(") and -vars- (") +
            std::to_string(vars.size()) +
            std::string(") should match. ")            
            );

    // Checking that there are variables in the model
    if (model != nullptr)
    {

        size_t K = model->get_agents_data_ncols();
        for (const auto & var: vars)
        {
            if ((var >= static_cast<int>(K)) | (var < 0))
                throw std::range_error(
                    std::string("The variable ") +
                    std::to_string(var) +
                    std::string(" is out of range.") +
                    std::string(" The agents only feature ") +
                    std::to_string(K) + 
                    std::string("variables (features).")
                );
        }
        
    }

    std::vector< epiworld_double > coefs_f;
    for (auto c: coefs)
        coefs_f.push_back(static_cast<epiworld_double>(c));

    ToolFun<TSeq> fun_ = [coefs_f,vars](
        Tool<TSeq>& tool,
        Agent<TSeq> * agent,
        VirusPtr<TSeq> virus,
        Model<TSeq> * model
        ) -> epiworld_double {

        size_t K = coefs_f.size();
        epiworld_double res = 0.0;

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd reduction(+:res)
        #endif
        for (size_t i = 0u; i < K; ++i)
            res += agent->operator[](vars.at(i)) * coefs_f.at(i);

        return 1.0/(1.0 + std::exp(-res));

    };

    return fun_;

}

template<typename TSeq>
inline Tool<TSeq>::Tool()
{
    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        sequence = -1;
    }
    else
    {
        sequence = nullptr;
    }

    set_name("Tool");
}

template<typename TSeq>
inline Tool<TSeq>::Tool(std::string name)
{
    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        sequence = -1;
    }
    else
    {
        sequence = nullptr;
    }

    set_name(name);
}

template<typename TSeq>
inline Tool<TSeq>::Tool(
    std::string name,
    epiworld_double prevalence,
    bool as_proportion
    )
{

    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        sequence = -1;
    }
    else
    {
        sequence = nullptr;
    }

    set_name(name);

    set_distribution(
        distribute_tool_randomly<TSeq>(prevalence, as_proportion)
    );
}

template<typename TSeq>
inline void Tool<TSeq>::set_sequence(TSeq d) {
    sequence = std::make_shared<TSeq>(d);
}

template<typename TSeq>
inline void Tool<TSeq>::set_sequence(std::shared_ptr<TSeq> d) {
    sequence = d;
}

template<>
inline void Tool<int>::set_sequence(int d) {
    sequence = d;
}

template<typename TSeq>
inline EPI_TYPENAME_TRAITS(TSeq, int) Tool<TSeq>::get_sequence() {
    return sequence;
}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (tool_functions->susceptibility_reduction)
        return tool_functions->susceptibility_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_CONTAGION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (tool_functions->transmission_reduction)
        return tool_functions->transmission_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_TRANSMISSION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (tool_functions->recovery_enhancer)
        return tool_functions->recovery_enhancer(*this, this->agent, v, model);

    return DEFAULT_TOOL_RECOVERY_ENHANCER;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_death_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (tool_functions->death_reduction)
        return tool_functions->death_reduction(*this, this->agent, v, model);

    return DEFAULT_TOOL_DEATH_REDUCTION;

}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction_fun(
    ToolFun<TSeq> fun
)
{
    tool_functions->susceptibility_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction_fun(
    ToolFun<TSeq> fun
)
{
    tool_functions->transmission_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer_fun(
    ToolFun<TSeq> fun
)
{
    tool_functions->recovery_enhancer = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction_fun(
    ToolFun<TSeq> fun
)
{
    tool_functions->death_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    tool_functions->susceptibility_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(susceptibility_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(epiworld_double * prob)
{
    
    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    tool_functions->transmission_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(transmission_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    tool_functions->recovery_enhancer = tmpfun;

}

// EPIWORLD_SET_LAMBDA(recovery_enhancer)
template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    tool_functions->death_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(death_reduction)

// #undef EPIWORLD_SET_LAMBDA
template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    tool_functions->susceptibility_reduction = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    tool_functions->transmission_reduction = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    tool_functions->recovery_enhancer = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    tool_functions->death_reduction = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_name(std::string name)
{
    tool_name = name;
}

template<typename TSeq>
inline std::string Tool<TSeq>::get_name() const {

    return tool_name;

}

template<typename TSeq>
inline Agent<TSeq> * Tool<TSeq>::get_agent()
{
    return this->agent;
}

template<typename TSeq>
inline void Tool<TSeq>::set_agent(Agent<TSeq> * p, size_t idx)
{
    agent        = p;
    pos_in_agent = static_cast<int>(idx);
}

template<typename TSeq>
inline int Tool<TSeq>::get_id() const {
    return id;
}


template<typename TSeq>
inline void Tool<TSeq>::set_id(int id)
{
    this->id = id;
}

template<typename TSeq>
inline void Tool<TSeq>::set_date(int d)
{
    this->date = d;
}

template<typename TSeq>
inline int Tool<TSeq>::get_date() const
{
    return date;
}

template<typename TSeq>
inline void Tool<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    state_init = init;
    state_post = end;
}

template<typename TSeq>
inline void Tool<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    queue_init = init;
    queue_post = end;
}

template<typename TSeq>
inline void Tool<TSeq>::get_state(
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
inline void Tool<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = queue_init;

    if (post != nullptr)
        *post = queue_post;

}

template<>
inline bool Tool<std::vector<int>>::operator==(
    const Tool<std::vector<int>> & other
    ) const
{
    
    if (sequence->size() != other.sequence->size())
        return false;

    for (size_t i = 0u; i < sequence->size(); ++i)
    {
        if (sequence->operator[](i) != other.sequence->operator[](i))
            return false;
    }

    if (tool_name != other.tool_name)
        return false;
    
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
inline bool Tool<TSeq>::operator==(const Tool<TSeq> & other) const
{
    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        if (sequence != other.sequence)
            return false;
    }
    else
    {
        if (*sequence != *other.sequence)
            return false;
    }


    if (tool_name != other.tool_name)
        return false;
    
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
inline void Tool<TSeq>::print() const
{

    printf_epiworld("Tool       : %s\n", tool_name.c_str());
    printf_epiworld("Id         : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());
    printf_epiworld("state_init : %i\n", static_cast<int>(state_init));
    printf_epiworld("state_post : %i\n", static_cast<int>(state_post));
    printf_epiworld("queue_init : %i\n", static_cast<int>(queue_init));
    printf_epiworld("queue_post : %i\n", static_cast<int>(queue_post));

}

template<typename TSeq>
inline void Tool<TSeq>::distribute(Model<TSeq> * model)
{

    if (tool_functions->dist)
    {

        tool_functions->dist(*this, model);

    }

}

template<typename TSeq>
inline void Tool<TSeq>::set_distribution(ToolToAgentFun<TSeq> fun)
{
    tool_functions->dist = fun;
}

#endif