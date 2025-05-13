#ifndef MEASLESQUARANTINE_HPP
#define MEASLESQUARANTINE_HPP

#if __cplusplus >= 202302L
    // C++23 or later
    #define GET_MODEL(model, output) \
        ModelMeaslesQuarantine<TSeq> * output = \
            dynamic_cast<ModelMeaslesQuarantine<TSeq> *>(model); \
        [[assume(output != nullptr)]]
#else
    // C++17 or C++20
    #define GET_MODEL(model, output) \
        ModelMeaslesQuarantine<TSeq> * output = \
            dynamic_cast<ModelMeaslesQuarantine<TSeq> *>(model); \
        assert(output != nullptr); // Use assert for runtime checks
#endif

#define LOCAL_UPDATE_FUN(name) \
    template<typename TSeq> \
    inline void ModelMeaslesQuarantine<TSeq>:: name \
    (epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m)

#define SAMPLE_FROM_PROBS(n, ans) \
    size_t ans; \
    epiworld_double p_total = m->runif(); \
    for (ans = 0u; ans < n; ++ans) \
    { \
        if (p_total < m->array_double_tmp[ans]) \
            break; \
        m->array_double_tmp[ans + 1] += m->array_double_tmp[ans]; \
    }

/**
 * @brief Template for a Measles model with quarantine
 * 
 * @param TSeq The type of the sequence to be used.
 * @details
 * This model can be described as a SEIHR model with isolation and quarantine.
 * The infectious state is divided into prodromal and rash phases. Furthermore,
 * the quarantine state includes exposed, susceptible, prodromal, and recovered
 * states.
 * 
 * The quarantine process is triggered any time that an agent with rash is
 * detected. The agent is then isolated and all agents who are unvaccinated are
 * quarantined. Isolated agents then may be moved out of the isolation in 
 * isolation_period days.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelMeaslesQuarantine: public Model<TSeq> {

private:

    /**
     * @brief The function that updates the model.
     */
    ///@{
    static void m_update_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_rash(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_q_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_q_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_q_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_q_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_hospitalized(Agent<TSeq> * p, Model<TSeq> * m);
    ///@}

    /**
     * @brief The function that updates the model.
     * 
     * This function is called at the end of each day.
     */
    static void m_update_model(Model<TSeq> * m);

public:

    static const epiworld_fast_uint SUSCEPTIBLE             = 0u;
    static const epiworld_fast_uint EXPOSED                 = 1u;
    static const epiworld_fast_uint PRODROMAL               = 2u;
    static const epiworld_fast_uint RASH                    = 3u;
    static const epiworld_fast_uint ISOLATED                = 4u;
    static const epiworld_fast_uint ISOLATED_RECOVERED      = 5u;
    static const epiworld_fast_uint DETECTED_HOSPITALIZED   = 6u;
    static const epiworld_fast_uint QUARANTINED_EXPOSED     = 7u;
    static const epiworld_fast_uint QUARANTINED_SUSCEPTIBLE = 8u;
    static const epiworld_fast_uint QUARANTINED_PRODROMAL   = 9u;
    static const epiworld_fast_uint QUARANTINED_RECOVERED   = 10u;
    static const epiworld_fast_uint HOSPITALIZED            = 11u;
    static const epiworld_fast_uint RECOVERED               = 12u;
    
    // Default constructor
    ModelMeaslesQuarantine() {};

    /**
     * @param n The number of agents in the system.
     * @param n_exposed The number of exposed agents in the system.
     * @param contact_rate The rate of contact between agents.
     * @param transmission_rate The rate of transmission of the virus.
     * @param vax_efficacy The efficacy of the vaccine.
     * @param vax_reduction_recovery_rate The reduction in recovery rate due to  the vaccine.
     * @param incubation_period The incubation period of the virus.
     * @param prodromal_period The prodromal period of the virus.
     * @param rash_period The rash period of the virus.
     * @param days_undetected The number of days the virus goes undetected.
     * @param hospitalization_rate The rate of hospitalization.
     * @param hospitalization_period The duration of hospitalization.
     * @param prop_vaccinated The proportion of vaccinated agents.
     * @param quarantine_period The number of days for quarantine.
     * @param quarantine_willingness The willingness to be quarantined.
     * @param isolation_period The number of days for isolation.
     */
    ///@{ 
    ModelMeaslesQuarantine(
        ModelMeaslesQuarantine<TSeq> & model,
        epiworld_fast_uint n,
        epiworld_fast_uint n_exposed,
        // Disease parameters
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double vax_efficacy,
        epiworld_double vax_reduction_recovery_rate,
        epiworld_double incubation_period,
        epiworld_double prodromal_period,
        epiworld_double rash_period,
        epiworld_double days_undetected,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double prop_vaccinated,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_fast_int isolation_period
    );

    ModelMeaslesQuarantine(
        epiworld_fast_uint n,
        epiworld_fast_uint n_exposed,
        // Disease parameters
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double vax_efficacy,
        epiworld_double vax_reduction_recovery_rate,
        epiworld_double incubation_period,
        epiworld_double prodromal_period,
        epiworld_double rash_period,
        epiworld_double days_undetected,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double prop_vaccinated,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_fast_int isolation_period
    );
    ///@}

    std::vector<Agent<TSeq> *> infectious; ///< Agents infectious for contact

    bool system_quarantine_triggered = false;

    std::vector< int > day_flagged; ///< Either detected or started quarantine
    std::vector< int > day_rash_onset; ///< Day of rash onset

    /**
     * @brief Quarantine agents that are in the system.
     * 
     * The flow should be:
     * - The function only runs if the quarantine status is active.
     * 
     * - Agents who are in quarantine, isolation, removed, or 
     *   hospitalized are ignored.
     * 
     * - Agents who are in the RASH state are isolated.
     * 
     * - Vaccinated agents are ignored.
     * 
     * - Susceptible, Exposed, and Prodromal agents are moved to the
     *   QUARANTINED_* state.
     * 
     * - At the end of the function, the quarantine status is set false.
     */
    void quarantine_agents();

    void reset();
    void update_infectious();

    Model<TSeq> * clone_ptr();

};

template<typename TSeq>
inline void ModelMeaslesQuarantine<TSeq>::quarantine_agents() {

    // Iterating through the new cases
    if (!system_quarantine_triggered)
        return;

    // Quarantine and isolation can be shut off if negative
    if (
        (this->par("Quarantine period") < 0) &&
        (this->par("Isolation period") < 0)
    )
        return;

    // Capturing the days that matter and the probability of success
    epiworld_double willingness = this->par("Quarantine willingness");

    // Iterating through the
    for (size_t i = 0u; i < this->size(); ++i) {

        auto agent_state = this->get_agent(i).get_state();

        // Already quarantined or isolated
        if (agent_state >= RASH)
            continue;

        // If the agent has a vaccine, then no need for quarantine
        if (this->get_agent(i).get_n_tools() != 0u)
            continue;

        // Quarantine will depend on the willingness of the agent
        // to be quarantined. If negative, then quarantine never happens.
        if (
            (this->par("Quarantine period") >= 0) &&
            (this->runif() < willingness)
        )
        {

            if (agent_state == SUSCEPTIBLE)
                this->get_agent(i).change_state(this, QUARANTINED_SUSCEPTIBLE);
            else if (agent_state == EXPOSED)
                this->get_agent(i).change_state(this, QUARANTINED_EXPOSED);
            else if (agent_state == PRODROMAL)
                this->get_agent(i).change_state(this, QUARANTINED_PRODROMAL);

            // And we add the day of quarantine
            this->day_flagged[i] = this->today();

        }

    }

    // Clearing the list of ids
    this->system_quarantine_triggered = false;

    return;

}


template<typename TSeq>
inline void ModelMeaslesQuarantine<TSeq>::m_update_model(Model<TSeq> * m) {
    
    GET_MODEL(m, model);
    model->quarantine_agents();
    model->events_run();
    model->update_infectious();
    return;

}

template<typename TSeq>
inline void ModelMeaslesQuarantine<TSeq>::reset() {
    
    Model<TSeq>::reset();

    this->system_quarantine_triggered = false;
        
    this->day_flagged.resize(this->size(), 0);
    std::fill(
        day_flagged.begin(),
        day_flagged.end(),
        0);

    this->day_rash_onset.resize(this->size(), 0);
    std::fill(
        day_rash_onset.begin(),
        day_rash_onset.end(),
        0);

    this->m_update_model(dynamic_cast<Model<TSeq>*>(this));
    return;
    
}

template<typename TSeq>
inline void ModelMeaslesQuarantine<TSeq>::update_infectious() {

    #ifdef EPI_DEBUG
    // All agents with state >= EXPOSED should have a virus
    for (auto & agent: this->get_agents())
    {
        auto s = agent.get_state();
        if (IN(s, {EXPOSED, PRODROMAL, RASH, ISOLATED, DETECTED_HOSPITALIZED, QUARANTINED_EXPOSED, QUARANTINED_PRODROMAL, HOSPITALIZED}))
        {
            if (agent.get_virus() == nullptr)
                throw std::logic_error("The agent has no virus.");
        }
    }
    #endif

    this->infectious.clear();
    int n_available = 0;
    for (auto & agent: this->get_agents())
    {
        const auto & s = agent.get_state();
        if (s == PRODROMAL)
            this->infectious.push_back(&agent);

        if (s < RASH)
            ++n_available;
        
    }

    // Assumes fixed contact rate throughout the simulation
    double p_contact = this->par("Contact rate")/
        static_cast< epiworld_double >(n_available);

    this->set_rand_binom(
        static_cast<int>(this->infectious.size()),
        p_contact > 1.0 ? 1.0 : p_contact
    );

}

template<typename TSeq>
inline Model<TSeq> * ModelMeaslesQuarantine<TSeq>::clone_ptr()
{
        
    ModelMeaslesQuarantine<TSeq> * ptr = new ModelMeaslesQuarantine<TSeq>(
        *dynamic_cast<const ModelMeaslesQuarantine<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

LOCAL_UPDATE_FUN(m_update_susceptible) {

    // How many contacts to draw
    int ndraw = m->rbinom();

    if (ndraw == 0)
        return;

    GET_MODEL(m, model);
    size_t n_infectious = model->infectious.size();

    if (n_infectious == 0)
        return;

    // Drawing from the set
    int nviruses_tmp = 0;
    int i = 0;
    while (i < ndraw)
    {
        // Picking the actual contacts
        int which = static_cast<int>(std::floor(n_infectious * m->runif()));

        /* There is a bug in which runif() returns 1.0. It is rare, but
            * we saw it here. See the Notes section in the C++ manual
            * https://en.cppreference.com/mwiki/index.php?title=cpp/numeric/random/uniform_real_distribution&oldid=133329
            * And the reported bug in GCC:
            * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
            * 
            */
        if (which == static_cast<int>(n_infectious))
            --which;

        epiworld::Agent<> & neighbor = *model->infectious[which];

        // Can't sample itself
        if (neighbor.get_id() == p->get_id())
            continue;

        // We successfully drew a contact, so we increment the counter
        i++;

        // No virus, the error!!
        if (neighbor.get_virus() == nullptr)
            throw std::logic_error("The neighbor has no virus.");

        // Only prodomal individuals can transmit
        if (neighbor.get_state() != model->PRODROMAL)
            throw std::logic_error(
                "The neighbor is not in the prodromal state. The state is: " +
                std::to_string(neighbor.get_state())
            );

        auto & v = neighbor.get_virus();

        #ifdef EPI_DEBUG
        if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
        #endif
            
        /* And it is a function of susceptibility_reduction as well */ 
        m->array_double_tmp[nviruses_tmp] =
            (1.0 - p->get_susceptibility_reduction(v, m)) * 
            v->get_prob_infecting(m) * 
            (1.0 - neighbor.get_transmission_reduction(v, m)) 
            ; 

        m->array_virus_tmp[nviruses_tmp++] = &(*v);
            
    }

    // No virus to compute
    if (nviruses_tmp == 0u)
        return;

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return;

    p->set_virus(*m->array_virus_tmp[which], m);

    return; 

};

LOCAL_UPDATE_FUN(m_update_exposed) {

    if (m->runif() < (1.0/p->get_virus()->get_incubation(m)))
        p->change_state(m, ModelMeaslesQuarantine<TSeq>::PRODROMAL);

    return;

};

LOCAL_UPDATE_FUN(m_update_prodromal) {
    
    if (m->runif() < (1.0/m->par("Prodromal period")))
    {

        GET_MODEL(m, model);
        model->day_rash_onset[p->get_id()] = m->today();
        p->change_state(m, ModelMeaslesQuarantine<TSeq>::RASH);

    }

    return;

};

LOCAL_UPDATE_FUN(m_update_rash) {


    GET_MODEL(m, model);
    
    #ifdef EPI_DEBUG
    if (model->day_flagged.size() <= p->get_id())
        throw std::logic_error(
            "The agent is not in the list of quarantined or isolated agents: " +
            std::to_string(p->get_id()) +
            " vs " +
            std::to_string(model->day_flagged.size()) +
            ". The model has " + std::to_string(model->size()) + " agents."
        );
    #endif

    // Checking if the agent will be detected or not
    // How many days since detected  
    bool detected = false;
    if (
        (m->par("Isolation period") >= 0) &&
        (m->runif() < 1.0/m->par("Days undetected"))
    )
    {
        model->system_quarantine_triggered = true;
        detected = true;

    }

    // Probability of Staying in the rash period vs becoming
    // hospitalized
    m->array_double_tmp[0] = 1.0/m->par("Rash period");
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    // Sampling from the probabilities
    SAMPLE_FROM_PROBS(2, which);

    // Recovers
    if (which == 2)
    {
        p->rm_agent_by_virus(
            m,
            detected ?
                ModelMeaslesQuarantine::ISOLATED_RECOVERED:
                ModelMeaslesQuarantine::RECOVERED
        );
    }
    else if (which == 1)
    {
        // If hospitalized, then the agent is removed from the system
        // effectively
        p->change_state(
            m,
            detected ?
                ModelMeaslesQuarantine::DETECTED_HOSPITALIZED :
                ModelMeaslesQuarantine::HOSPITALIZED
            );
    }
    else if (which != 0)
    {
        throw std::logic_error("The roulette returned an unexpected value.");
    } else if ((which == 0u) && detected)
    {
        // If the agent is not hospitalized, then it is moved to
        // isolation.
        p->change_state(m, ModelMeaslesQuarantine::ISOLATED);
    }
    
};

LOCAL_UPDATE_FUN(m_update_isolated) {

    GET_MODEL(m, model);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_rash_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    // Probability of staying in the rash period vs becoming
    // hospitalized
    m->array_double_tmp[0] = 1.0/m->par("Rash period");
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    // Sampling from the probabilities
    SAMPLE_FROM_PROBS(2, which);

    // Recovers
    if (which == 2u)
    {
        if (unisolate)
        {
            p->rm_agent_by_virus(
                m,
                ModelMeaslesQuarantine::RECOVERED
            );
        }
        else
            p->rm_agent_by_virus(
                m, ModelMeaslesQuarantine::ISOLATED_RECOVERED
            );
    }

    // If hospitalized, then the agent is removed from the system
    else if (which == 1u)
    {
        p->change_state(m, ModelMeaslesQuarantine::HOSPITALIZED);
    }
    // If neither hospitalized nor recovered, then the agent is
    // still under isolation, unless the quarantine period is over.
    else if ((which == 0u) && unisolate)
    {
        p->change_state(m, ModelMeaslesQuarantine::RASH);
    }

}

LOCAL_UPDATE_FUN(m_update_isolated_recovered) {

    GET_MODEL(m, model);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_rash_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    if (unisolate)
        p->change_state(m, ModelMeaslesQuarantine::RECOVERED);

}

LOCAL_UPDATE_FUN(m_update_q_exposed) {

    // How many days since quarantine started
    GET_MODEL(m, model);
    int days_since =
        m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    // Will develop prodromal symptoms?
    if (m->runif() < (1.0/p->get_virus()->get_incubation(m)))
    {
        // If the quarantine period is over, then they are moved to
        // the prodromal period. Otherwise, they are moved to the
        // quarantined prodromal period.
        if (unquarantine)
            p-> change_state(
                m,
                ModelMeaslesQuarantine::PRODROMAL
            );
        else
            p->change_state(
                m,
                ModelMeaslesQuarantine::QUARANTINED_PRODROMAL
            );

    }
    else if (unquarantine)
    {
        p->change_state(
            m,
            ModelMeaslesQuarantine::EXPOSED
        );
    }

}

LOCAL_UPDATE_FUN(m_update_q_susceptible) {

    GET_MODEL(m, model);
    int days_since =
        m->today() - model->day_flagged[p->get_id()];
    
    if (days_since >= m->par("Quarantine period"))
        p->change_state(m, ModelMeaslesQuarantine::SUSCEPTIBLE);

}

LOCAL_UPDATE_FUN(m_update_q_prodromal) {

    GET_MODEL(m, model);

    // Otherwise, these are moved to the prodromal period, if
    // the quanrantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;
    
    // Develops rash?
    if (m->runif() < (1.0/m->par("Prodromal period")))
    {
        model->day_rash_onset[p->get_id()] = m->today();
        p->change_state(m, ModelMeaslesQuarantine::ISOLATED);
    }
    else
    {
        
        if (unquarantine)
            p->change_state(m, ModelMeaslesQuarantine::PRODROMAL);

    }

}

LOCAL_UPDATE_FUN(m_update_q_recovered) {

    GET_MODEL(m, model);
    int days_since = m->today() - model->day_flagged[p->get_id()];
    
    if (days_since >= m->par("Quarantine period"))
        p->change_state(m, ModelMeaslesQuarantine::RECOVERED);

}

LOCAL_UPDATE_FUN(m_update_hospitalized) {

    // The agent is removed from the system
    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_agent_by_virus(m, ModelMeaslesQuarantine::RECOVERED);

    return;

}


template<typename TSeq>
inline ModelMeaslesQuarantine<TSeq>::ModelMeaslesQuarantine(
    ModelMeaslesQuarantine<TSeq> & model,
    epiworld_fast_uint n,
    epiworld_fast_uint n_exposed,
    // Disease parameters
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double vax_efficacy,
    epiworld_double vax_reduction_recovery_rate,
    epiworld_double incubation_period,
    epiworld_double prodromal_period,
    epiworld_double rash_period,
    epiworld_double days_undetected,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double prop_vaccinated,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_fast_int isolation_period
) {

    model.add_state("Susceptible", this->m_update_susceptible);
    model.add_state("Exposed", this->m_update_exposed);
    model.add_state("Prodromal", this->m_update_prodromal);
    model.add_state("Rash", this->m_update_rash);
    model.add_state("Isolated", this->m_update_isolated);
    model.add_state(
        "Isolated Recovered", this->m_update_isolated_recovered
    );
    model.add_state("Detected Hospitalized", this->m_update_hospitalized);
    model.add_state(
        "Quarantined Exposed", this->m_update_q_exposed
    );

    model.add_state(
        "Quarantined Susceptible", this->m_update_q_susceptible
    );

    model.add_state(
        "Quarantined Prodromal", this->m_update_q_prodromal
    );

    model.add_state(
        "Quarantined Recovered", this->m_update_q_recovered
    );

    model.add_state("Hospitalized", this->m_update_hospitalized);

    model.add_state("Recovered");

    // Adding the model parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(incubation_period, "Incubation period");
    model.add_param(prodromal_period, "Prodromal period");
    model.add_param(rash_period, "Rash period");
    model.add_param(days_undetected, "Days undetected");
    model.add_param(quarantine_period, "Quarantine period");
    model.add_param(
        quarantine_willingness, "Quarantine willingness"
    );
    model.add_param(isolation_period, "Isolation period");
    model.add_param(hospitalization_rate, "Hospitalization rate");
    model.add_param(hospitalization_period, "Hospitalization period");
    model.add_param(prop_vaccinated, "Vaccination rate");
    model.add_param(vax_efficacy, "Vax efficacy");
    model.add_param(vax_reduction_recovery_rate, "Vax improved recovery");

    // Designing the disease
    Virus<> measles("Measles");
    measles.set_state(EXPOSED, RECOVERED);
    measles.set_prob_infecting(&model("Transmission rate"));
    measles.set_prob_recovery(&model("Rash period"));
    measles.set_incubation(&model("Incubation period"));
    measles.set_distribution(
        distribute_virus_randomly(n_exposed, false)
    );

    model.add_virus(measles);

    // Designing the vaccine
    Tool<> vaccine("Vaccine");
    vaccine.set_susceptibility_reduction(&model("Vax efficacy"));
    vaccine.set_recovery_enhancer(&model("Vax improved recovery"));
    vaccine.set_distribution(
        distribute_tool_randomly(prop_vaccinated, true)
    );

    model.add_tool(vaccine);

    // Global actions
    model.add_globalevent(this->m_update_model, "Update model");
    model.queuing_off();

    // Setting the population
    model.agents_empty_graph(n);

    return;

}

template<typename TSeq>
inline ModelMeaslesQuarantine<TSeq>::ModelMeaslesQuarantine(
    epiworld_fast_uint n,
    epiworld_fast_uint n_exposed,
    // Disease parameters
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double vax_efficacy,
    epiworld_double vax_reduction_recovery_rate,
    epiworld_double incubation_period,
    epiworld_double prodromal_period,
    epiworld_double rash_period,
    epiworld_double days_undetected,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double prop_vaccinated,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_fast_int isolation_period

) {

    ModelMeaslesQuarantine(
        *this,
        n,
        n_exposed,
        contact_rate,
        transmission_rate,
        vax_efficacy,
        vax_reduction_recovery_rate,
        incubation_period,
        prodromal_period,
        rash_period,
        days_undetected,
        hospitalization_rate,
        hospitalization_period,
        prop_vaccinated,
        quarantine_period,
        quarantine_willingness,
        isolation_period
    );

    return;

}

#undef SAMPLE_FROM_PROBS
#undef LOCAL_UPDATE_FUN
#undef GET_MODEL
#endif