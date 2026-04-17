#ifndef EPIWORLD_MODELS_SEIRNETWORKQUARANTINE_HPP
#define EPIWORLD_MODELS_SEIRNETWORKQUARANTINE_HPP

#include "../model-bones.hpp"

/**
 * @file seirnetworkquarantine.hpp
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 *        with network-based contacts, quarantine, and contact tracing
 */

/**
 * @brief SEIR model with network-based contacts, quarantine, and contact tracing
 *
 * This class implements a Susceptible-Exposed-Infected-Removed (SEIR) epidemiological
 * model that uses actual network neighbors (e.g., from a Stochastic Block Model) for
 * contacts instead of a mixing matrix. This enables the queueing system for efficient
 * iteration since only agents in contact with infected agents need updates.
 *
 * Features:
 * - Network-based contacts (via adjacency list, e.g., from agents_sbm())
 * - Quarantine measures for exposed contacts
 * - Isolation policies for detected infected individuals
 * - Contact tracing with configurable success rates
 * - Hospitalization of severe cases
 * - Individual willingness to comply with public health measures
 * - Queueing system support for computational efficiency
 *
 * The model supports 10 distinct states (identical to SEIRMixingQuarantine):
 * - Susceptible: Individuals who can become infected
 * - Exposed: Infected but not yet infectious (incubation period)
 * - Infected: Infectious individuals in the community
 * - Isolated: Detected infected individuals in self-isolation
 * - Detected Hospitalized: Hospitalized individuals who were contact-traced
 * - Quarantined Susceptible: Susceptible individuals in quarantine due to contact tracing
 * - Quarantined Exposed: Exposed individuals in quarantine due to contact tracing
 * - Isolated Recovered: Recovered individuals still in isolation
 * - Hospitalized: Individuals requiring hospital care
 * - Recovered: Individuals who have recovered and gained immunity
 *
 * @tparam TSeq Type for genetic sequences (default: EPI_DEFAULT_TSEQ)
 * @ingroup connected_models
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRNetworkQuarantine : public Model<TSeq>
{
private:

    // Update functions
    static void _update_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_infected(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_isolated(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_quarantine_suscep(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_quarantine_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_hospitalized(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_isolated_recovered(Agent<TSeq> * p, Model<TSeq> * m);

    // Data about the quarantine process
    std::vector< bool > quarantine_willingness; ///< Indicator
    std::vector< bool > isolation_willingness; ///< Indicator for isolation willingness
    std::vector< size_t > agent_quarantine_triggered; ///< Whether the quarantine process has started
    std::vector< int > day_flagged; ///< Either detected or started quarantine
    std::vector< int > day_onset; ///< Day of onset of the disease

    static void _quarantine_process(Model<TSeq> * m);

public:

    static const int SUSCEPTIBLE             = 0;
    static const int EXPOSED                 = 1;
    static const int INFECTED                = 2;
    static const int ISOLATED                = 3;
    static const int DETECTED_HOSPITALIZED   = 4;
    static const int QUARANTINED_SUSCEPTIBLE = 5;
    static const int QUARANTINED_EXPOSED     = 6;
    static const int ISOLATED_RECOVERED      = 7;
    static const int HOSPITALIZED            = 8;
    static const int RECOVERED               = 9;

    static const size_t QUARANTINE_PROCESS_INACTIVE = 0u;
    static const size_t QUARANTINE_PROCESS_ACTIVE   = 1u;
    static const size_t QUARANTINE_PROCESS_DONE     = 2u;

    ModelSEIRNetworkQuarantine() = delete;

    /**
     * @brief Constructs a ModelSEIRNetworkQuarantine object.
     *
     * @param vname The name of the virus.
     * @param n The number of agents in the model.
     * @param prevalence The initial prevalence of the disease.
     * @param transmission_rate The transmission rate of the disease.
     * @param avg_incubation_days The average incubation period.
     * @param recovery_rate The recovery rate of the disease.
     * @param hospitalization_rate The rate at which infected individuals are hospitalized.
     * @param hospitalization_period The average duration of hospitalization in days.
     * @param days_undetected The average number of days an infected individual remains undetected.
     * @param quarantine_period The duration of quarantine in days.
     * @param quarantine_willingness The proportion willing to comply with quarantine.
     * @param isolation_willingness The proportion willing to self-isolate.
     * @param isolation_period The duration of isolation in days.
     * @param contact_tracing_success_rate Probability of successfully tracing a contact (default: 1.0).
     * @param contact_tracing_days_prior Days prior to detection for contact tracing (default: 4).
     */
    ModelSEIRNetworkQuarantine(
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double days_undetected,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_double isolation_willingness,
        epiworld_fast_int isolation_period,
        epiworld_double contact_tracing_success_rate = 1.0,
        epiworld_fast_uint contact_tracing_days_prior = 4u
    );

    /**
     * @brief Reset the model to initial state
     */
    void reset() override;

    /**
     * @brief Create a clone of this model
     * @return Pointer to a new model instance with the same configuration
     */
    std::unique_ptr< Model<TSeq> > clone_ptr() override;

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with two elements:
     * - [0]: The proportion of initially infected individuals who start in the exposed state.
     * - [1]: The proportion of initially non-infected individuals who have recovered (immune).
     * @param queue_ Optional vector for queuing specifications (default: empty).
     */
    ModelSEIRNetworkQuarantine<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    ) override;

    /**
     * @brief Get the quarantine trigger status for all agents
     * @return Vector indicating quarantine process status for each agent
     */
    std::vector< size_t > get_agent_quarantine_triggered() const
    {
        return agent_quarantine_triggered;
    };

    /**
     * @brief Get the quarantine willingness for all agents
     * @return Vector of boolean values indicating each agent's willingness to quarantine
     */
    std::vector< bool > get_quarantine_willingness() const
    {
        return quarantine_willingness;
    };

    /**
     * @brief Get the isolation willingness for all agents
     * @return Vector of boolean values indicating each agent's willingness to self-isolate
     */
    std::vector< bool > get_isolation_willingness() const
    {
        return isolation_willingness;
    };

};

// -----------------------------------------------------------------------
// Reset
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::reset()
{
    Model<TSeq>::reset();

    // Setting up the quarantine parameters
    quarantine_willingness.resize(this->size(), false);
    isolation_willingness.resize(this->size(), false);
    for (size_t idx = 0; idx < quarantine_willingness.size(); ++idx)
    {
        quarantine_willingness[idx] =
            this->runif() < this->par("Quarantine willingness");
        isolation_willingness[idx] =
            this->runif() < this->par("Isolation willingness");
    }

    agent_quarantine_triggered.assign(this->size(), 0u);
    day_flagged.assign(this->size(), 0);
    day_onset.assign(this->size(), 0);

    return;
}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSEIRNetworkQuarantine<TSeq>::clone_ptr()
{
    return std::make_unique<ModelSEIRNetworkQuarantine<TSeq>>(*this);
}

// -----------------------------------------------------------------------
// Susceptible: iterate over network neighbors (like default_update_susceptible)
// and record contacts for tracing
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_susceptible(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    size_t nviruses_tmp = 0u;
    for (auto & neighbor : p->get_neighbors(*m))
    {
        auto & v = neighbor->get_virus();
        if (v == nullptr)
            continue;

        // Only infectious agents can transmit
        if (neighbor->get_state() != ModelSEIRNetworkQuarantine<TSeq>::INFECTED)
            continue;

        // Record contact for tracing: infected neighbor -> susceptible agent
        m->get_contact_tracing().add_contact(
            neighbor->get_id(),
            p->get_id(),
            static_cast<size_t>(m->today())
        );

        #ifdef EPI_DEBUG
        if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
            throw std::logic_error(
                "Trying to add an extra element to a temporal array outside of the range."
            );
        #endif

        m->array_double_tmp[nviruses_tmp] =
            (1.0 - p->get_susceptibility_reduction(v, *m)) *
            v->get_prob_infecting(m) *
            (1.0 - neighbor->get_transmission_reduction(v, *m));

        m->array_virus_tmp[nviruses_tmp++] = &(*v);
    }

    if (nviruses_tmp == 0u)
        return;

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return;

    p->set_virus(*m,
        *m->array_virus_tmp[which],
        ModelSEIRNetworkQuarantine<TSeq>::EXPOSED
    );

    return;
};

// -----------------------------------------------------------------------
// Exposed: incubation -> infected transition
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto & v = p->get_virus();

    if (m->runif() < 1.0/(v->get_incubation(m)))
    {
        p->change_state(*m, ModelSEIRNetworkQuarantine<TSeq>::INFECTED);

        auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);
        model->day_onset[p->get_id()] = m->today();
    }

    return;
};

// -----------------------------------------------------------------------
// Infected: detection, isolation, hospitalization, recovery
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_infected(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    // Sampling whether the agent is detected or not.
    // If Days undetected < 0, detection is disabled (never detected).
    // If Days undetected == 0, the agent is always detected.
    epiworld_double days_undetected = m->par("Days undetected");
    bool detected = (days_undetected < 0.0) ?
        false : ((days_undetected == 0.0) ?
            true : (m->runif() < 1.0 / days_undetected));

    // If detected, trigger the quarantine process
    if (detected)
    {
        model->agent_quarantine_triggered[p->get_id()] =
            ModelSEIRNetworkQuarantine<TSeq>::QUARANTINE_PROCESS_ACTIVE;
    }

    // Checking if the agent is willing to isolate individually
    bool isolation_detected = (m->par("Isolation period") >= 0) &&
        detected &&
        (model->isolation_willingness[p->get_id()]);

    // Recording the date of detection
    if (isolation_detected)
        model->day_flagged[p->get_id()] = m->today();

    // Computing probabilities for state change
    auto & v = p->get_virus();
    m->array_double_tmp[0] = 1.0 - (1.0 - v->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(v, *m));
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    auto which = m->sample_from_probs(2);

    if (which == 0) // Recovers
    {
        if (isolation_detected)
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::ISOLATED_RECOVERED
            );
        }
        else
        {
            p->rm_virus(*m,
                ModelSEIRNetworkQuarantine<TSeq>::RECOVERED
            );
        }
        return;
    }
    else if (which == 1) // Hospitalized
    {
        if (detected)
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::DETECTED_HOSPITALIZED
            );
        }
        else
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::HOSPITALIZED
            );
        }
    }
    else if ((which == 2) && isolation_detected) // Nothing, but detected
    {
        p->change_state(*m,
            ModelSEIRNetworkQuarantine<TSeq>::ISOLATED
        );
    }

    return;
};

// -----------------------------------------------------------------------
// Isolated: recovery, hospitalization, or release from isolation
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_isolated(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true : false;

    // Sampling from the probabilities of recovery
    m->array_double_tmp[0] = 1.0 -
        (1.0 - p->get_virus()->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(p->get_virus(), *m));

    // And hospitalization
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    auto which = m->sample_from_probs(2);

    // Recovers
    if (which == 0)
    {
        if (unisolate)
        {
            p->rm_virus(*m,
                ModelSEIRNetworkQuarantine<TSeq>::RECOVERED
            );
        }
        else
            p->rm_virus(*m,
                ModelSEIRNetworkQuarantine<TSeq>::ISOLATED_RECOVERED
            );
    }
    else if (which == 1)
    {
        if (unisolate)
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::HOSPITALIZED
            );
        }
        else
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::DETECTED_HOSPITALIZED
            );
        }
    }
    else if ((which == 2) && unisolate)
    {
        p->change_state(*m,
            ModelSEIRNetworkQuarantine<TSeq>::INFECTED
        );
    }

};

// -----------------------------------------------------------------------
// Quarantined Susceptible: release when quarantine period is over
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_quarantine_suscep(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true : false;

    if (unquarantine)
    {
        p->change_state(*m,
            ModelSEIRNetworkQuarantine<TSeq>::SUSCEPTIBLE
        );
    }

};

// -----------------------------------------------------------------------
// Quarantined Exposed: incubation or release
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_quarantine_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true : false;

    if (m->runif() < 1.0/(p->get_virus()->get_incubation(m)))
    {
        // Recording the day of onset
        model->day_onset[p->get_id()] = m->today();

        if (unquarantine)
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::INFECTED
            );
        }
        else
        {
            p->change_state(*m,
                ModelSEIRNetworkQuarantine<TSeq>::ISOLATED
            );
        }
    }
    else if (unquarantine)
    {
        p->change_state(*m,
            ModelSEIRNetworkQuarantine<TSeq>::EXPOSED
        );
    }

};

// -----------------------------------------------------------------------
// Isolated Recovered: release when isolation period ends
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_isolated_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true : false;

    if (unisolate)
    {
        p->change_state(*m,
            ModelSEIRNetworkQuarantine<TSeq>::RECOVERED
        );
    }

};

// -----------------------------------------------------------------------
// Hospitalized: recovery after hospitalization period
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_update_hospitalized(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_virus(*m, ModelSEIRNetworkQuarantine<TSeq>::RECOVERED);

};

// -----------------------------------------------------------------------
// Quarantine process: trace contacts and quarantine them
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRNetworkQuarantine<TSeq>::_quarantine_process(
    Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRNetworkQuarantine<TSeq>, TSeq>(m);

    for (size_t agent_i = 0u; agent_i < m->size(); ++agent_i)
    {

        if (
            model->agent_quarantine_triggered[agent_i] !=
            ModelSEIRNetworkQuarantine<TSeq>::QUARANTINE_PROCESS_ACTIVE
        )
            continue;

        if (m->par("Quarantine period") < 0)
        {
            model->agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;
            continue;
        }

        size_t n_contacts = m->get_contact_tracing().get_n_contacts(agent_i);
        if (n_contacts >= EPI_MAX_TRACKING)
            n_contacts = EPI_MAX_TRACKING;

        auto success_rate = m->par("Contact tracing success rate");
        auto days_prior = m->par("Contact tracing days prior");
        for (size_t contact_i = 0u; contact_i < n_contacts; ++contact_i)
        {
            // Checking if we will detect the contact
            if (m->runif() > success_rate)
                continue;

            auto [contact_id, contact_date] = m->get_contact_tracing().get_contact(
                agent_i, contact_i
            );

            // Skip contacts outside the tracing window
            if ((static_cast<double>(m->today()) -
                static_cast<double>(contact_date)) > days_prior)
                continue;

            auto & agent = m->get_agent(contact_id);

            if (agent.get_state() > INFECTED)
                continue;

            // Agents with some tool won't be quarantined
            if (agent.get_n_tools() != 0u)
                continue;

            if (model->quarantine_willingness[contact_id])
            {
                switch (agent.get_state())
                {
                    case SUSCEPTIBLE:
                        agent.change_state(*m, QUARANTINED_SUSCEPTIBLE);
                        model->day_flagged[contact_id] = m->today();
                        break;
                    case EXPOSED:
                        agent.change_state(*m, QUARANTINED_EXPOSED);
                        model->day_flagged[contact_id] = m->today();
                        break;
                    case INFECTED:
                        if (model->isolation_willingness[contact_id])
                        {
                            agent.change_state(*m, ISOLATED);
                            model->day_flagged[contact_id] = m->today();
                        }
                        break;
                    default:
                        throw std::logic_error(
                            "The agent is not in a state that can be quarantined."
                        );
                }
            }
        }

        // Setting the quarantine process off
        model->agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;
    }

    return;
}

// -----------------------------------------------------------------------
// Constructor (delegating)
// -----------------------------------------------------------------------
template<typename TSeq>
inline ModelSEIRNetworkQuarantine<TSeq>::ModelSEIRNetworkQuarantine(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double days_undetected,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_double isolation_willingness,
    epiworld_fast_int isolation_period,
    epiworld_double contact_tracing_success_rate,
    epiworld_fast_uint contact_tracing_days_prior
    )
{

    // Setting up parameters
    this->add_param(transmission_rate, "Prob. Transmission");
    this->add_param(recovery_rate, "Prob. Recovery");
    this->add_param(avg_incubation_days, "Avg. Incubation days");
    this->add_param(hospitalization_rate, "Hospitalization rate");
    this->add_param(hospitalization_period, "Hospitalization period");
    this->add_param(days_undetected, "Days undetected");
    this->add_param(quarantine_period, "Quarantine period");
    this->add_param(
        quarantine_willingness, "Quarantine willingness"
    );
    this->add_param(
        isolation_willingness, "Isolation willingness"
    );
    this->add_param(isolation_period, "Isolation period");
    this->add_param(
        contact_tracing_success_rate, "Contact tracing success rate"
    );
    this->add_param(
        contact_tracing_days_prior, "Contact tracing days prior"
    );

    // States
    this->add_state("Susceptible", _update_susceptible);
    this->add_state("Exposed", _update_exposed);
    this->add_state("Infected", _update_infected);
    this->add_state("Isolated", _update_isolated);
    this->add_state("Detected Hospitalized", _update_hospitalized);
    this->add_state("Quarantined Susceptible", _update_quarantine_suscep);
    this->add_state("Quarantined Exposed", _update_quarantine_exposed);
    this->add_state("Isolated Recovered", _update_isolated_recovered);
    this->add_state("Hospitalized", _update_hospitalized);
    this->add_state("Recovered");

    // Global function (quarantine process runs before state updates)
    this->add_globalevent(_quarantine_process, "Quarantine process");

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(EXPOSED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Prob. Transmission");
    virus.set_prob_recovery("Prob. Recovery");
    virus.set_incubation("Avg. Incubation days");

    this->add_virus(virus);

    // Enable contact tracing for quarantine process
    this->contact_tracing_on(EPI_MAX_TRACKING);

    this->set_name("SEIR with Network and Quarantine");

    return;

}

template<typename TSeq>
inline ModelSEIRNetworkQuarantine<TSeq> & ModelSEIRNetworkQuarantine<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /* queue_ */
)
{

    Model<TSeq>::initial_states_fun =
        create_init_function_seir<TSeq>(proportions_)
        ;

    return *this;

}

#endif
