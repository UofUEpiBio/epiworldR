#ifndef EPIWORLD_VIRUS_MEAT_HPP
#define EPIWORLD_VIRUS_MEAT_HPP

/**
 * @brief Factory function of VirusFun base on logit
 * 
 * @tparam TSeq 
 * @param vars Vector indicating the position of the variables to use.
 * @param coefs Vector of coefficients.
 * @return VirusFun<TSeq> 
 */
template<typename TSeq>
inline VirusFun<TSeq> virus_fun_logit(
    std::vector< int > vars,
    std::vector< double > coefs,
    Model<TSeq> * model,
    bool logit = true
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

    VirusFun<TSeq> fun_infect = [coefs_f,vars](
        Agent<TSeq> * agent,
        Virus<TSeq> &,
        Model<TSeq> *
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

    return fun_infect;

}

template<typename TSeq>
inline Virus<TSeq>::Virus(
    std::string name
    ) {
    set_name(name);
}

template<typename TSeq>
inline Virus<TSeq>::Virus(
    std::string name,
    epiworld_double prevalence,
    bool prevalence_as_proportion
    ) {
    set_name(name);
    set_distribution(
        distribute_virus_randomly<TSeq>(
            prevalence,
            prevalence_as_proportion
        )
    );
}

template<typename TSeq>
inline void Virus<TSeq>::mutate(
    Model<TSeq> * model
) {

    if (mutation_fun)
        if (mutation_fun(agent, *this, model))
            model->get_db().record_virus(*this);

    return;
    
}

template<typename TSeq>
inline void Virus<TSeq>::set_mutation(
    MutFun<TSeq> fun
) {
    mutation_fun = MutFun<TSeq>(fun);
}

template<typename TSeq>
inline std::shared_ptr<TSeq> Virus<TSeq>::get_sequence()
{

    return baseline_sequence;

}

template<typename TSeq>
inline void Virus<TSeq>::set_sequence(TSeq sequence)
{

    baseline_sequence = std::make_shared<TSeq>(sequence);
    return;

}

template<typename TSeq>
inline Agent<TSeq> * Virus<TSeq>::get_agent()
{

    return agent;

}

template<typename TSeq>
inline void Virus<TSeq>::set_agent(Agent<TSeq> * p)
{
    agent = p;
}

template<typename TSeq>
inline void Virus<TSeq>::set_id(int idx)
{

    id = idx;
    return;

}

template<typename TSeq>
inline int Virus<TSeq>::get_id() const
{
    
    return id;

}

template<typename TSeq>
inline void Virus<TSeq>::set_date(int d) 
{

    date = d;
    return;

}

template<typename TSeq>
inline int Virus<TSeq>::get_date() const
{
    
    return date;
    
}

template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_infecting(
    Model<TSeq> * model
)
{

    if (probability_of_infecting_fun)
        return probability_of_infecting_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_INFECTION;

}



template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_recovery(
    Model<TSeq> * model
)
{

    if (probability_of_recovery_fun)
        return probability_of_recovery_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_RECOVERY;

}



template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_death(
    Model<TSeq> * model
)
{

    if (probability_of_death_fun)
        return probability_of_death_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_DEATH;

}

template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_incubation(
    Model<TSeq> * model
)
{

    if (incubation_fun)
        return incubation_fun(agent, *this, model);
        
    return EPI_DEFAULT_INCUBATION_DAYS;

}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting_fun(VirusFun<TSeq> fun)
{
    probability_of_infecting_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery_fun(VirusFun<TSeq> fun)
{
    probability_of_recovery_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death_fun(VirusFun<TSeq> fun)
{
    probability_of_death_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_incubation_fun(VirusFun<TSeq> fun)
{
    incubation_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting(const epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_infecting_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery(const epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_recovery_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death(const epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_death_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_incubation(const epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    incubation_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_infecting_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_recovery_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_death_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_incubation(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    incubation_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_post_recovery(PostRecoveryFun<TSeq> fun)
{
    if (post_recovery_fun)
    {
        printf_epiworld(
            "Warning: a PostRecoveryFun is alreay in place (overwriting)."
            );
    }

    post_recovery_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::post_recovery(
    Model<TSeq> * model
)
{

    if (post_recovery_fun)
        post_recovery_fun(agent, *this, model);    

    return;
        
}

template<typename TSeq>
inline void Virus<TSeq>::set_post_immunity(
    epiworld_double prob
)
{

    if (post_recovery_fun)
    {

        std::string msg =
            std::string(
                "You cannot set post immunity when a post_recovery "
                ) +
            std::string(
                "function is already in place. Redesign the post_recovery function."
                );

        throw std::logic_error(msg);
        
    }

    // To make sure that we keep registering the virus
    ToolPtr<TSeq> __no_reinfect = std::make_shared<Tool<TSeq>>(
        "Immunity (" + *virus_name + ")"
    );

    __no_reinfect->set_susceptibility_reduction(prob);
    __no_reinfect->set_death_reduction(0.0);
    __no_reinfect->set_transmission_reduction(0.0);
    __no_reinfect->set_recovery_enhancer(0.0);

    PostRecoveryFun<TSeq> tmpfun = 
        [__no_reinfect](
            Agent<TSeq> * p, Virus<TSeq> &, Model<TSeq> * m
            )
        {
            
            // Have we registered the tool?
            if (__no_reinfect->get_id() == -99)
                m->get_db().record_tool(*__no_reinfect);

            p->add_tool(__no_reinfect, m);

            return;

        };

    post_recovery_fun = tmpfun;

}

template<typename TSeq>
inline void Virus<TSeq>::set_post_immunity(
    epiworld_double * prob
)
{

    if (post_recovery_fun)
    {

        std::string msg =
            std::string(
                "You cannot set post immunity when a post_recovery "
                ) +
            std::string(
                "function is already in place. Redesign the post_recovery function."
                );

        throw std::logic_error(msg);

    }

    // To make sure that we keep registering the virus
    ToolPtr<TSeq> __no_reinfect = std::make_shared<Tool<TSeq>>(
        "Immunity (" + *virus_name + ")"
    );

    __no_reinfect->set_susceptibility_reduction(prob);
    __no_reinfect->set_death_reduction(0.0);
    __no_reinfect->set_transmission_reduction(0.0);
    __no_reinfect->set_recovery_enhancer(0.0);

    PostRecoveryFun<TSeq> tmpfun = 
        [__no_reinfect](Agent<TSeq> * p, Virus<TSeq> &, Model<TSeq> * m)
        {

            // Have we registered the tool?
            if (__no_reinfect->get_id() == -99)
                m->get_db().record_tool(*__no_reinfect);

            p->add_tool(__no_reinfect, m);

            return;

        };

    post_recovery_fun = tmpfun;

}

template<typename TSeq>
inline void Virus<TSeq>::set_name(std::string name)
{

    if (name == "")
        virus_name = nullptr;
    else
        virus_name = std::make_shared<std::string>(name);

}

template<typename TSeq>
inline std::string Virus<TSeq>::get_name() const
{

    if (virus_name)
        return *virus_name;
    
    return "unknown virus";

}

template<typename TSeq>
inline std::vector< epiworld_double > & Virus<TSeq>::get_data() {
    return data;
}

template<typename TSeq>
inline void Virus<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end,
    epiworld_fast_int removed
)
{
    state_init    = init;
    state_post    = end;
    state_removed = removed;
}

template<typename TSeq>
inline void Virus<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end,
    epiworld_fast_int removed
)
{

    queue_init    = init;
    queue_post     = end;
    queue_removed = removed;

}

template<typename TSeq>
inline void Virus<TSeq>::get_state(
    epiworld_fast_int * init,
    epiworld_fast_int * end,
    epiworld_fast_int * removed
)
{

    if (init != nullptr)
        *init = state_init;

    if (end != nullptr)
        *end = state_post;

    if (removed != nullptr)
        *removed = state_removed;

}

template<typename TSeq>
inline void Virus<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * end,
    epiworld_fast_int * removed
)
{

    if (init != nullptr)
        *init = queue_init;

    if (end != nullptr)
        *end = queue_post;

    if (removed != nullptr)
        *removed = queue_removed;
        
}

template<>
inline bool Virus<std::vector<int>>::operator==(
    const Virus<std::vector<int>> & other
    ) const
{
    

    EPI_DEBUG_FAIL_AT_TRUE(
        baseline_sequence->size() != other.baseline_sequence->size(),
        "Virus:: baseline_sequence don't match"
        )

    for (size_t i = 0u; i < baseline_sequence->size(); ++i)
    {

        EPI_DEBUG_FAIL_AT_TRUE(
            baseline_sequence->operator[](i) != other.baseline_sequence->operator[](i),
            "Virus:: baseline_sequence[i] don't match"
            )

    }

    EPI_DEBUG_FAIL_AT_TRUE(
        virus_name != other.virus_name,
        "Virus:: virus_name don't match"
        )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        state_init != other.state_init,
        "Virus:: state_init don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        state_post != other.state_post,
        "Virus:: state_post don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        state_removed != other.state_removed,
        "Virus:: state_removed don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_init != other.queue_init,
        "Virus:: queue_init don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_post != other.queue_post,
        "Virus:: queue_post don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_removed != other.queue_removed,
        "Virus:: queue_removed don't match"
        )

    return true;

}

template<typename TSeq>
inline bool Virus<TSeq>::operator==(const Virus<TSeq> & other) const
{
    
    EPI_DEBUG_FAIL_AT_TRUE(
        *baseline_sequence != *other.baseline_sequence,
        "Virus:: baseline_sequence don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        virus_name != other.virus_name,
        "Virus:: virus_name don't match"
    )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        state_init != other.state_init,
        "Virus:: state_init don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        state_post != other.state_post,
        "Virus:: state_post don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        state_removed != other.state_removed,
        "Virus:: state_removed don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_init != other.queue_init,
        "Virus:: queue_init don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_post != other.queue_post,
        "Virus:: queue_post don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_removed != other.queue_removed,
        "Virus:: queue_removed don't match"
    )

    return true;

}

template<typename TSeq>
inline void Virus<TSeq>::print() const
{

    printf_epiworld("Virus         : %s\n", virus_name->c_str());
    printf_epiworld("Id            : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());
    printf_epiworld("state_init    : %i\n", static_cast<int>(state_init));
    printf_epiworld("state_post    : %i\n", static_cast<int>(state_post));
    printf_epiworld("state_removed : %i\n", static_cast<int>(state_removed));
    printf_epiworld("queue_init    : %i\n", static_cast<int>(queue_init));
    printf_epiworld("queue_post    : %i\n", static_cast<int>(queue_post));
    printf_epiworld("queue_removed : %i\n", static_cast<int>(queue_removed));

}

template<typename TSeq>
inline void Virus<TSeq>::distribute(Model<TSeq> * model)
{

    if (dist_fun)
    {

        dist_fun(*this, model);

    }

}

template<typename TSeq>
inline void Virus<TSeq>::set_distribution(VirusToAgentFun<TSeq> fun)
{
    dist_fun = fun;
}

#endif
