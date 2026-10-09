
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
        Tool<TSeq>&,
        Agent<TSeq> * agent,
        VirusPtr<TSeq> &,
        Model<TSeq> * model
        ) -> epiworld_double {

        size_t K = coefs_f.size();
        epiworld_double res = 0.0;

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd reduction(+:res)
        #endif
        for (size_t i = 0u; i < K; ++i)
            res += agent->operator()(vars.at(i), *model) * coefs_f.at(i);

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
    VirusPtr<TSeq> & v,
    Model<TSeq> * model
)
{

    if (susceptibility_reduction)
        return susceptibility_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_CONTAGION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> * model
)
{

    if (transmission_reduction)
        return transmission_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_TRANSMISSION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> & v,
    Model<TSeq> * model
)
{

    if (recovery_enhancer)
        return recovery_enhancer(*this, this->agent, v, model);

    return DEFAULT_TOOL_RECOVERY_ENHANCER;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_death_reduction(
    VirusPtr<TSeq> & v,
    Model<TSeq> * model
)
{

    if (death_reduction)
        return death_reduction(*this, this->agent, v, model);

    return DEFAULT_TOOL_DEATH_REDUCTION;

}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction_fun(
    ToolFun<TSeq> fun
)
{
    susceptibility_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction_fun(
    ToolFun<TSeq> fun
)
{
    transmission_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer_fun(
    ToolFun<TSeq> fun
)
{
    recovery_enhancer = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction_fun(
    ToolFun<TSeq> fun
)
{
    death_reduction = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun<TSeq> tmpfun =
        [param_ref](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq>* model)
        {
            return (*param_ref)(*model);
        };

    susceptibility_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(susceptibility_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));
    
    ToolFun<TSeq> tmpfun =
        [param_ref](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq>* model)
        {
            return (*param_ref)(*model);
        };

    transmission_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(transmission_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun<TSeq> tmpfun =
        [param_ref](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq>* model)
        {
            return (*param_ref)(*model);
        };

    recovery_enhancer = tmpfun;

}

// EPIWORLD_SET_LAMBDA(recovery_enhancer)
template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun<TSeq> tmpfun =
        [param_ref](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq>* model)
        {
            return (*param_ref)(*model);
        };

    death_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(death_reduction)

// #undef EPIWORLD_SET_LAMBDA
template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq> *)
        {
            return prob;
        };

    susceptibility_reduction = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq> *)
        {
            return prob;
        };

    transmission_reduction = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq> *)
        {
            return prob;
        };

    recovery_enhancer = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>&, Model<TSeq> *)
        {
            return prob;
        };

    death_reduction = tmpfun;

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
inline uint64_t Tool<TSeq>::target_bit(int lineage_id)
{

    if ((lineage_id < 0) || (lineage_id >= 63))
        throw std::range_error(
            std::string("The virus lineage id ") +
            std::to_string(lineage_id) +
            std::string(" cannot be targeted. Only lineages 0 to 62 can be ") +
            std::string("targeted by tools.")
        );

    return uint64_t(1) << lineage_id;

}

template<typename TSeq>
inline void Tool<TSeq>::add_target(int lineage_id)
{

    uint64_t bit = target_bit(lineage_id);

    // The first target replaces the default (every virus)
    if (target_mask == ~uint64_t(0))
        target_mask = 0u;

    target_mask |= bit;

}

template<typename TSeq>
inline void Tool<TSeq>::add_target(const Virus<TSeq> & v)
{

    if (v.get_lineage_id() < 0)
        throw std::logic_error(
            std::string("The virus \"") + v.get_name() +
            std::string("\" has no lineage id. Add it to the model with ") +
            std::string("Model::add_virus() before targeting it.")
        );

    add_target(v.get_lineage_id());

}

template<typename TSeq>
inline void Tool<TSeq>::set_targets(const std::vector< int > & lineage_ids)
{

    // Validate every id before touching the current targets
    uint64_t mask = lineage_ids.empty() ? ~uint64_t(0) : 0u;
    for (auto id : lineage_ids)
        mask |= target_bit(id);

    target_mask = mask;

}

template<typename TSeq>
inline std::vector< int > Tool<TSeq>::get_targets() const
{

    std::vector< int > res;
    if (target_mask == ~uint64_t(0))
        return res;

    for (int i = 0; i < 63; ++i)
        if (target_mask & (uint64_t(1) << i))
            res.push_back(i);

    return res;

}

template<typename TSeq>
inline void Tool<TSeq>::clear_targets()
{
    target_mask = ~uint64_t(0);
}

template<typename TSeq>
inline bool Tool<TSeq>::targets(const Virus<TSeq> & v) const
{
    return (target_mask & v.lineage_bit) != 0u;
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

    if (target_mask != other.target_mask)
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

    if (target_mask != other.target_mask)
        return false;

    return true;

}


template<typename TSeq>
inline void Tool<TSeq>::print() const
{

    printf_epiworld("Tool       : %s\n", this->get_name().c_str());
    printf_epiworld("Id         : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());
    printf_epiworld("state_init : %i\n", static_cast<int>(state_init));
    printf_epiworld("state_post : %i\n", static_cast<int>(state_post));
    printf_epiworld("queue_init : %i\n", static_cast<int>(queue_init));
    printf_epiworld("queue_post : %i\n", static_cast<int>(queue_post));

    if (target_mask != ~uint64_t(0))
    {
        std::string tgts;
        for (auto i : get_targets())
            tgts += (tgts.empty() ? "" : ", ") + std::to_string(i);
        printf_epiworld("targets    : %s\n", tgts.c_str());
    }

}

template<typename TSeq>
inline void Tool<TSeq>::distribute(Model<TSeq> * model)
{

    if (dist)
    {

        dist(*this, model);

    }

}

template<typename TSeq>
inline void Tool<TSeq>::set_distribution(ToolToAgentFun<TSeq> fun)
{
    dist = fun;
}

template<typename TSeq>
inline std::unique_ptr<Tool<TSeq>> Tool<TSeq>::clone_ptr() const
{
    auto cloned = std::make_unique<Tool<TSeq>>(*this);
    return cloned;
}

#endif
