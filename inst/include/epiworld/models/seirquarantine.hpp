#ifndef EPIWORLD_MODELS_SEIRQUARANTINE_HPP
#define EPIWORLD_MODELS_SEIRQUARANTINE_HPP

#include "../model-bones.hpp"

/**
 * @file seirquarantine.hpp
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 *        with quarantine, isolation, and contact tracing
 */

/**
 * @brief SEIR model with quarantine, isolation, and contact tracing
 *
 * @details An SEIR model whose contacts come from its layers, of any kind
 * (see `Model::add_layer()`): networks, mixing by group, or both. Features:
 *
 * - Detection of infected agents, who isolate if willing.
 * - Contact tracing from the contact logs of the model's layers: at reset the
 *   model turns tracing on for each layer that does not have it yet (see
 *   `ContactLayer::enable_tracing()`), and a detected agent's contacts in every
 *   log are traced.
 * - Quarantine of the traced contacts who are willing to comply.
 * - Hospitalization of severe cases.
 * - Individual willingness to quarantine and to isolate (drawn once per agent
 *   and run).
 *
 * Only Infected agents transmit (`Model::set_infectious_states()`). Isolated,
 * quarantined and hospitalized agents are out of the mixing layers altogether
 * (`Model::set_excluded_states()`): they are not met, and their group's
 * contacts go to the others.
 *
 * The model has 10 states:
 *
 * - Susceptible: Individuals who can become infected
 * - Exposed: Infected but not yet infectious (latent period)
 * - Infected: Infectious individuals in the community
 * - Isolated: Detected infected individuals in self-isolation
 * - Detected Hospitalized: Hospitalized individuals who were detected
 * - Quarantined Susceptible: Susceptible individuals in quarantine due to contact tracing
 * - Quarantined Exposed: Exposed individuals in quarantine due to contact tracing
 * - Isolated Recovered: Recovered individuals still in isolation
 * - Hospitalized: Individuals requiring hospital care
 * - Recovered: Individuals who have recovered and gained immunity
 *
 * **Implementation details:**
 * <a href="../impl/quarantine-isolation-and-contact-tracing.md">Quarantine, Isolation, and Contact Tracing</a>
 *
 * @tparam TSeq Type for genetic sequences (default: EPI_DEFAULT_TSEQ)
 * @ingroup basic_compartmental
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRQuarantine : public Model<TSeq>
{
private:

    // Update functions
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

    ModelSEIRQuarantine() = delete;

    /**
     * @brief Constructs a ModelSEIRQuarantine object.
     *
     * @param vname The name of the virus.
     * @param n The number of agents in the model.
     * @param prevalence The initial prevalence of the disease.
     * @param transmission_prob Probability that a contact with an infected
     * agent infects ("Transmission probability").
     * @param latent_period Mean days in Exposed ("Latent period").
     * @param recovery_rate Daily probability of recovering ("Recovery rate").
     * @param hospitalization_rate Daily probability of being hospitalized
     * while infected ("Hospitalization rate").
     * @param hospitalization_period Mean days in hospital ("Hospitalization
     * period").
     * @param detection_delay Mean days infected before detection ("Detection
     * delay"; 0 detects at once, a negative value never).
     * @param quarantine_period Days in quarantine ("Quarantine period"; a
     * negative value turns quarantine off).
     * @param quarantine_willingness Probability that an agent complies with
     * quarantine ("Quarantine willingness").
     * @param isolation_willingness Probability that an agent isolates when
     * detected ("Isolation willingness").
     * @param isolation_period Days in isolation, from onset ("Isolation
     * period"; a negative value turns isolation off).
     * @param contact_tracing_prob Probability that a contact is traced
     * ("Contact tracing probability", default 1).
     * @param contact_tracing_window Days before detection whose contacts are
     * traced ("Contact tracing window", default 4).
     */
    ModelSEIRQuarantine(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_prob,
        epiworld_double latent_period,
        epiworld_double recovery_rate,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double detection_delay,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_double isolation_willingness,
        epiworld_fast_int isolation_period,
        epiworld_double contact_tracing_prob = 1.0,
        epiworld_fast_uint contact_tracing_window = 4u
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
    ModelSEIRQuarantine<TSeq> & initial_states(
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
inline void ModelSEIRQuarantine<TSeq>::reset()
{

    // The quarantine process reads the contacts of every layer
    for (size_t i = 0u; i < this->get_n_layers(); ++i)
        if (!this->get_layer(i).is_tracing())
            this->get_layer(i).enable_tracing();

    Model<TSeq>::reset();

    // Setting up the quarantine parameters
    quarantine_willingness.resize(this->size(), false);
    isolation_willingness.resize(this->size(), false);
    for (size_t idx = 0; idx < quarantine_willingness.size(); ++idx)
    {
        quarantine_willingness[idx] =
            this->runif() < EPI_PAR(this, "Quarantine willingness");
        isolation_willingness[idx] =
            this->runif() < EPI_PAR(this, "Isolation willingness");
    }

    agent_quarantine_triggered.assign(this->size(), 0u);
    day_flagged.assign(this->size(), 0);
    day_onset.assign(this->size(), 0);

    return;
}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSEIRQuarantine<TSeq>::clone_ptr()
{
    return std::make_unique<ModelSEIRQuarantine<TSeq>>(*this);
}

// -----------------------------------------------------------------------
// Exposed: incubation -> infected transition
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto & v = p->get_virus();

    if (m->runif() < 1.0/(v->get_incubation(m)))
    {
        p->change_state(*m, ModelSEIRQuarantine<TSeq>::INFECTED);

        auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);
        model->day_onset[p->get_id()] = m->today();
    }

    return;
};

// -----------------------------------------------------------------------
// Infected: detection, isolation, hospitalization, recovery
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_infected(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    // Sampling whether the agent is detected or not.
    // A negative detection delay disables detection (never detected); 0
    // detects at once.
    epiworld_double detection_delay = EPI_PAR(m, "Detection delay");
    bool detected = (detection_delay < 0.0) ?
        false : ((detection_delay == 0.0) ?
            true : (m->runif() < 1.0 / detection_delay));

    // If detected, trigger the quarantine process
    if (detected)
    {
        model->agent_quarantine_triggered[p->get_id()] =
            ModelSEIRQuarantine<TSeq>::QUARANTINE_PROCESS_ACTIVE;
    }

    // Checking if the agent is willing to isolate individually
    bool isolation_detected = (EPI_PAR(m, "Isolation period") >= 0) &&
        detected &&
        (model->isolation_willingness[p->get_id()]);

    // Recording the date of detection
    if (isolation_detected)
        model->day_flagged[p->get_id()] = m->today();

    // Computing probabilities for state change
    auto & v = p->get_virus();
    m->array_double_tmp[0] = 1.0 - (1.0 - v->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(v, *m));
    m->array_double_tmp[1] = EPI_PAR(m, "Hospitalization rate");

    auto which = m->sample_from_probs(2);

    if (which == 0) // Recovers
    {
        if (isolation_detected)
        {
            p->change_state(*m,
                ModelSEIRQuarantine<TSeq>::ISOLATED_RECOVERED
            );
        }
        else
        {
            p->rm_virus(*m,
                ModelSEIRQuarantine<TSeq>::RECOVERED
            );
        }
        return;
    }
    else if (which == 1) // Hospitalized
    {
        m->record_hospitalization(*p);
        if (detected)
        {
            p->change_state(*m, DETECTED_HOSPITALIZED);
        }
        else
        {
            p->change_state(*m, HOSPITALIZED);
        }
    }
    else if ((which == 2) && isolation_detected) // Nothing, but detected
    {
        p->change_state(*m, ISOLATED);
    }

    return;
};

// -----------------------------------------------------------------------
// Isolated: recovery, hospitalization, or release from isolation
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_isolated(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (EPI_PAR(m, "Isolation period") <= days_since) ?
        true : false;

    // Sampling from the probabilities of recovery
    m->array_double_tmp[0] = 1.0 -
        (1.0 - p->get_virus()->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(p->get_virus(), *m));

    // And hospitalization
    m->array_double_tmp[1] = EPI_PAR(m, "Hospitalization rate");

    auto which = m->sample_from_probs(2);

    // Recovers
    if (which == 0)
    {
        if (unisolate)
        {
            p->rm_virus(*m,
                ModelSEIRQuarantine<TSeq>::RECOVERED
            );
        }
        else
            p->rm_virus(*m,
                ModelSEIRQuarantine<TSeq>::ISOLATED_RECOVERED
            );
    }
    else if (which == 1)
    {
        m->record_hospitalization(*p);
        if (unisolate)
        {
            p->change_state(*m, HOSPITALIZED);
        }
        else
        {
            p->change_state(*m, DETECTED_HOSPITALIZED);
        }
    }
    else if ((which == 2) && unisolate)
    {
        p->change_state(*m, INFECTED);
    }

};

// -----------------------------------------------------------------------
// Quarantined Susceptible: release when quarantine period is over
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_quarantine_suscep(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (EPI_PAR(m, "Quarantine period") <= days_since) ?
        true : false;

    if (unquarantine)
    {
        p->change_state(*m,
            ModelSEIRQuarantine<TSeq>::SUSCEPTIBLE
        );
    }

};

// -----------------------------------------------------------------------
// Quarantined Exposed: incubation or release
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_quarantine_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (EPI_PAR(m, "Quarantine period") <= days_since) ?
        true : false;

    if (m->runif() < 1.0/(p->get_virus()->get_incubation(m)))
    {
        // Recording the day of onset
        model->day_onset[p->get_id()] = m->today();

        if (unquarantine)
        {
            p->change_state(*m,
                ModelSEIRQuarantine<TSeq>::INFECTED
            );
        }
        else
        {
            p->change_state(*m,
                ModelSEIRQuarantine<TSeq>::ISOLATED
            );
        }
    }
    else if (unquarantine)
    {
        p->change_state(*m,
            ModelSEIRQuarantine<TSeq>::EXPOSED
        );
    }

};

// -----------------------------------------------------------------------
// Isolated Recovered: release when isolation period ends
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_isolated_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (EPI_PAR(m, "Isolation period") <= days_since) ?
        true : false;

    if (unisolate)
    {
        p->change_state(*m,
            ModelSEIRQuarantine<TSeq>::RECOVERED
        );
    }

};

// -----------------------------------------------------------------------
// Hospitalized: recovery after hospitalization period
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_update_hospitalized(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    if (m->runif() < 1.0/EPI_PAR(m, "Hospitalization period"))
        p->rm_virus(*m, ModelSEIRQuarantine<TSeq>::RECOVERED);

};

// -----------------------------------------------------------------------
// Quarantine process: trace contacts and quarantine them
// -----------------------------------------------------------------------
template<typename TSeq>
inline void ModelSEIRQuarantine<TSeq>::_quarantine_process(
    Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRQuarantine<TSeq>, TSeq>(m);

    for (size_t agent_i = 0u; agent_i < m->size(); ++agent_i)
    {

        if (
            model->agent_quarantine_triggered[agent_i] !=
            ModelSEIRQuarantine<TSeq>::QUARANTINE_PROCESS_ACTIVE
        )
            continue;

        if (EPI_PAR(m, "Quarantine period") < 0)
        {
            model->agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;
            continue;
        }

        auto success_rate = EPI_PAR(m, "Contact tracing probability");
        auto days_prior = EPI_PAR(m, "Contact tracing window");

        // The contacts in the log of every traced layer, in the order the
        // layers were added (a log keeps each agent's last contacts only)
        for (size_t layer_i = 0u; layer_i < m->get_n_layers(); ++layer_i)
        {

            auto & layer = m->get_layer(layer_i);
            if (!layer.is_tracing())
                continue;

            auto & ct = layer.get_tracing();
            const size_t n_contacts = std::min(
                ct.get_n_contacts(agent_i), ct.get_max_contacts()
            );

            for (size_t contact_i = 0u; contact_i < n_contacts; ++contact_i)
            {
                // Checking if we will detect the contact
                if (m->runif() > success_rate)
                    continue;

                auto [contact_id, contact_date] = ct.get_contact(
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

        }

        // Setting the quarantine process off
        model->agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;
    }

    return;
}

// -----------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------
template<typename TSeq>
inline ModelSEIRQuarantine<TSeq>::ModelSEIRQuarantine(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double latent_period,
    epiworld_double recovery_rate,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double detection_delay,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_double isolation_willingness,
    epiworld_fast_int isolation_period,
    epiworld_double contact_tracing_prob,
    epiworld_fast_uint contact_tracing_window
    )
{

    this->set_population_size(n);

    // Setting up parameters
    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(latent_period, "Latent period");
    this->add_param(hospitalization_rate, "Hospitalization rate");
    this->add_param(hospitalization_period, "Hospitalization period");
    this->add_param(detection_delay, "Detection delay");
    this->add_param(quarantine_period, "Quarantine period");
    this->add_param(
        quarantine_willingness, "Quarantine willingness"
    );
    this->add_param(
        isolation_willingness, "Isolation willingness"
    );
    this->add_param(isolation_period, "Isolation period");
    this->add_param(contact_tracing_prob, "Contact tracing probability");
    this->add_param(contact_tracing_window, "Contact tracing window");

    // States
    this->add_state("Susceptible", default_update_susceptible<TSeq>);
    this->add_state("Exposed", _update_exposed);
    this->add_state("Infected", _update_infected);
    this->add_state("Isolated", _update_isolated);
    this->add_state("Detected Hospitalized", _update_hospitalized);
    this->add_state("Quarantined Susceptible", _update_quarantine_suscep);
    this->add_state("Quarantined Exposed", _update_quarantine_exposed);
    this->add_state("Isolated Recovered", _update_isolated_recovered);
    this->add_state("Hospitalized", _update_hospitalized);
    this->add_state("Recovered");

    // Only infected agents transmit; isolated, quarantined and hospitalized
    // ones are out of the mixing
    this->set_infectious_states({INFECTED});
    this->set_excluded_states({
        ISOLATED, DETECTED_HOSPITALIZED, QUARANTINED_SUSCEPTIBLE,
        QUARANTINED_EXPOSED, ISOLATED_RECOVERED, HOSPITALIZED
    });

    // Global function (quarantine process runs before state updates)
    this->add_globalevent(_quarantine_process, "Quarantine process");

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(EXPOSED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Transmission probability");
    virus.set_prob_recovery("Recovery rate");
    virus.set_incubation("Latent period");

    this->add_virus(virus);

    this->set_name("SEIR with Quarantine");

    return;

}

template<typename TSeq>
inline ModelSEIRQuarantine<TSeq> & ModelSEIRQuarantine<TSeq>::initial_states(
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
