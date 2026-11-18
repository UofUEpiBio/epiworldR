
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

inline void DiseaseTargets::add(std::string disease)
{
    if (std::find(names.begin(), names.end(), disease) == names.end())
        names.push_back(std::move(disease));
}

inline bool DiseaseTargets::has(const Virus & v) const
{
    return (mask & v.disease_bit) != 0u;
}

inline void Tool::add_target(std::string disease)
{
    targets_.add(std::move(disease));
}

inline void Tool::set_targets(std::vector< std::string > diseases)
{
    targets_ = DiseaseTargets();
    for (auto & d : diseases)
        targets_.add(std::move(d));
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
    




    if (!(targets_ == other.targets_))
        return false;

    return true;

}


inline void Tool::print() const
{

    printf_epiworld("Tool       : %s\n", this->get_name().c_str());
    printf_epiworld("Id         : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());

    if (!targets_.names.empty())
    {
        std::string tgts;
        for (const auto & n : targets_.names)
            tgts += (tgts.empty() ? "" : ", ") + n;
        printf_epiworld("targets    : %s\n", tgts.c_str());
    }

}

inline std::unique_ptr<Tool> Tool::clone_ptr() const
{
    auto cloned = std::make_unique<Tool>(*this);
    return cloned;
}

#endif
