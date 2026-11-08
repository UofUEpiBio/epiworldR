
#ifndef EPIWORLD_TOOLS_MEAT_HPP
#define EPIWORLD_TOOLS_MEAT_HPP

/**
 * @brief Factory function of ToolFun base on logit
 * 
 * @param vars Vector indicating the position of the variables to use.
 * @param coefs Vector of coefficients.
 * @return ToolFun 
 */
inline ToolFun tool_fun_logit(
    std::vector< int > vars,
    std::vector< double > coefs,
    Model * model
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

    ToolFun fun_ = [coefs_f,vars](
        Tool&,
        Agent * agent,
        VirusPtr &,
        Model * model
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

inline Tool::Tool()
{
    set_name("Tool");
}

inline Tool::Tool(std::string name)
{
    set_name(name);
}

inline epiworld_double Tool::get_susceptibility_reduction(
    VirusPtr & v,
    Model * model
)
{

    if (susceptibility_reduction)
        return susceptibility_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_CONTAGION_REDUCTION;

}

inline epiworld_double Tool::get_transmission_reduction(
    VirusPtr & v,
    Model * model
)
{

    if (transmission_reduction)
        return transmission_reduction(
            *this, this->agent, v, model
        );

    return DEFAULT_TOOL_TRANSMISSION_REDUCTION;

}

inline epiworld_double Tool::get_recovery_enhancer(
    VirusPtr & v,
    Model * model
)
{

    if (recovery_enhancer)
        return recovery_enhancer(*this, this->agent, v, model);

    return DEFAULT_TOOL_RECOVERY_ENHANCER;

}

inline epiworld_double Tool::get_death_reduction(
    VirusPtr & v,
    Model * model
)
{

    if (death_reduction)
        return death_reduction(*this, this->agent, v, model);

    return DEFAULT_TOOL_DEATH_REDUCTION;

}

inline void Tool::set_susceptibility_reduction_fun(
    ToolFun fun
)
{
    susceptibility_reduction = fun;
}

inline void Tool::set_transmission_reduction_fun(
    ToolFun fun
)
{
    transmission_reduction = fun;
}

inline void Tool::set_recovery_enhancer_fun(
    ToolFun fun
)
{
    recovery_enhancer = fun;
}

inline void Tool::set_death_reduction_fun(
    ToolFun fun
)
{
    death_reduction = fun;
}

inline void Tool::set_susceptibility_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun tmpfun =
        [param_ref](Tool &, Agent *, VirusPtr&, Model* model)
        {
            return (*param_ref)(*model);
        };

    susceptibility_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(susceptibility_reduction)
inline void Tool::set_transmission_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));
    
    ToolFun tmpfun =
        [param_ref](Tool &, Agent *, VirusPtr&, Model* model)
        {
            return (*param_ref)(*model);
        };

    transmission_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(transmission_reduction)
inline void Tool::set_recovery_enhancer(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun tmpfun =
        [param_ref](Tool &, Agent *, VirusPtr&, Model* model)
        {
            return (*param_ref)(*model);
        };

    recovery_enhancer = tmpfun;

}

// EPIWORLD_SET_LAMBDA(recovery_enhancer)
inline void Tool::set_death_reduction(std::string param)
{

    auto param_ref = std::make_shared< const ParamRef >(std::move(param));

    ToolFun tmpfun =
        [param_ref](Tool &, Agent *, VirusPtr&, Model* model)
        {
            return (*param_ref)(*model);
        };

    death_reduction = tmpfun;

}

// EPIWORLD_SET_LAMBDA(death_reduction)

// #undef EPIWORLD_SET_LAMBDA
inline void Tool::set_susceptibility_reduction(
    epiworld_double prob
)
{

    ToolFun tmpfun = 
        [prob](Tool &, Agent *, VirusPtr&, Model *)
        {
            return prob;
        };

    susceptibility_reduction = tmpfun;

}

inline void Tool::set_transmission_reduction(
    epiworld_double prob
)
{

    ToolFun tmpfun = 
        [prob](Tool &, Agent *, VirusPtr&, Model *)
        {
            return prob;
        };

    transmission_reduction = tmpfun;

}

inline void Tool::set_recovery_enhancer(
    epiworld_double prob
)
{

    ToolFun tmpfun = 
        [prob](Tool &, Agent *, VirusPtr&, Model *)
        {
            return prob;
        };

    recovery_enhancer = tmpfun;

}

inline void Tool::set_death_reduction(
    epiworld_double prob
)
{

    ToolFun tmpfun = 
        [prob](Tool &, Agent *, VirusPtr&, Model *)
        {
            return prob;
        };

    death_reduction = tmpfun;

}

inline void Tool::set_name(std::string name)
{
    tool_name = name;
}

inline std::string Tool::get_name() const {

    return tool_name;

}

inline uint64_t Tool::target_bit(int lineage_id)
{

    if ((lineage_id < 0) || (lineage_id >= 63))
        throw std::range_error(
            std::string("The virus lineage id ") +
            std::to_string(lineage_id) +
            std::string(" cannot be targeted. Only lineages 0 to 62 can be ") +
            std::string("targeted.")
        );

    return uint64_t(1) << lineage_id;

}

inline void Tool::add_target(int lineage_id)
{

    uint64_t bit = target_bit(lineage_id);

    // The first target replaces the default (every virus)
    if (target_mask == ~uint64_t(0))
        target_mask = 0u;

    target_mask |= bit;

}

inline void Tool::add_target(const Virus & v)
{

    if (v.get_lineage_id() < 0)
        throw std::logic_error(
            std::string("The virus \"") + v.get_name() +
            std::string("\" has no lineage id. Add it to the model with ") +
            std::string("Model::add_virus() before targeting it.")
        );

    add_target(v.get_lineage_id());

}

inline void Tool::set_targets(const std::vector< int > & lineage_ids)
{

    // Validate every id before touching the current targets
    uint64_t mask = lineage_ids.empty() ? ~uint64_t(0) : 0u;
    for (auto id : lineage_ids)
        mask |= target_bit(id);

    target_mask = mask;

}

inline std::vector< int > Tool::get_targets() const
{

    std::vector< int > res;
    if (target_mask == ~uint64_t(0))
        return res;

    for (int i = 0; i < 63; ++i)
        if (target_mask & (uint64_t(1) << i))
            res.push_back(i);

    return res;

}

inline void Tool::clear_targets()
{
    target_mask = ~uint64_t(0);
}

inline bool Tool::targets(const Virus & v) const
{
    return (target_mask & v.lineage_bit) != 0u;
}

inline Agent * Tool::get_agent()
{
    return this->agent;
}

inline void Tool::set_agent(Agent * p, size_t idx)
{
    agent        = p;
    pos_in_agent = static_cast<int>(idx);
}

inline int Tool::get_id() const {
    return id;
}


inline void Tool::set_id(int id)
{
    this->id = id;
}

inline void Tool::set_date(int d)
{
    this->date = d;
}

inline int Tool::get_date() const
{
    return date;
}

inline bool Tool::operator==(const Tool & other) const
{

    if (tool_name != other.tool_name)
        return false;
    




    if (target_mask != other.target_mask)
        return false;

    return true;

}


inline void Tool::print() const
{

    printf_epiworld("Tool       : %s\n", this->get_name().c_str());
    printf_epiworld("Id         : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());

    if (target_mask != ~uint64_t(0))
    {
        std::string tgts;
        for (auto i : get_targets())
            tgts += (tgts.empty() ? "" : ", ") + std::to_string(i);
        printf_epiworld("targets    : %s\n", tgts.c_str());
    }

}

inline std::unique_ptr<Tool> Tool::clone_ptr() const
{
    auto cloned = std::make_unique<Tool>(*this);
    return cloned;
}

#endif
