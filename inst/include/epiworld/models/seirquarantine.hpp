#ifndef EPIWORLD_MODELS_SEIRQUARANTINE_HPP
#define EPIWORLD_MODELS_SEIRQUARANTINE_HPP

#include "../model-bones.hpp"

/**
 * @file seirquarantine.hpp
 * @brief A Susceptible-Exposed-Infected-Removed (SEIR) model
 *        with quarantine, isolation, and contact tracing
 */

/**
 * @brief SEIR model with quarantine, isolation, and contact tracing
 *
 * ![Model Diagram](../assets/img/seirquarantine.png)
 *
 * @details An SEIR model whose contacts come from its layers, of any kind
 * (see `Model::add_layer()`): networks, mixing by group, or both. Each agent
 * has a stage of the disease and a status:
 *
 * - Stages (the virus's natural history): Susceptible, Exposed (latent),
 *   Infected (infectious), Severe (needs hospital care), Recovered. Infected
 *   agents recover or become severe (competing daily rates); severe ones
 *   recover at rate 1 / "Hospitalization period".
 * - Statuses (the model's): Free, Quarantined, Isolated, Hospitalized. The
 *   last three take the agent out of every layer.
 *
 * Every day, the quarantine process (a global event):
 *
 * 1. detects each free infected agent with probability 1 / "Detection delay",
 *    from the day after its onset (a negative delay disables detection, 0
 *    detects on the first day possible); a detected
 *    agent willing to isolate is isolated ("Isolated", if "Isolation period"
 *    is not negative);
 * 2. traces the contacts of the detected agents in every layer's log (the
 *    model turns tracing on at reset), within "Contact tracing window" days,
 *    each found with probability "Contact tracing probability". A free
 *    contact without tools who is willing to quarantine is quarantined if
 *    susceptible or exposed ("Quarantined"), or isolated if infected and
 *    willing to isolate ("Isolated").
 *
 * Declared rules:
 *
 * - Quarantine ends "Quarantine period" days after it starts.
 * - Isolation ends "Isolation period" days after the onset (entering
 *   Infected), whatever the stage then.
 * - An agent quarantined at the start of the day of its onset is isolated
 *   that day ("Onset in quarantine").
 * - Entering Severe admits the agent ("Admitted", recorded as a
 *   hospitalization); leaving it discharges the agent ("Discharged").
 *
 * **Implementation details:**
 * <a href="../impl/quarantine-isolation-and-contact-tracing.md">Quarantine, Isolation, and Contact Tracing</a>
 *
 * @ingroup basic_compartmental
 */
class ModelSEIRQuarantine : public Model
{
private:

    // Data about the quarantine process
    std::vector< bool > quarantine_willingness; ///< Indicator
    std::vector< bool > isolation_willingness; ///< Indicator for isolation willingness
    std::vector< size_t > agent_quarantine_triggered; ///< Whether the quarantine process has started

    static void _quarantine_process(Model * m);

public:

    // Stages
    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int EXPOSED     = 1;
    static constexpr int INFECTED    = 2;
    static constexpr int SEVERE      = 3;
    static constexpr int RECOVERED   = 4;

    // Statuses
    static constexpr int FREE         = 0;
    static constexpr int QUARANTINED  = 1;
    static constexpr int ISOLATED     = 2;
    static constexpr int HOSPITALIZED = 3;

    static const size_t QUARANTINE_PROCESS_INACTIVE = 0u;
    static const size_t QUARANTINE_PROCESS_ACTIVE   = 1u;
    static const size_t QUARANTINE_PROCESS_DONE     = 2u;

    ModelSEIRQuarantine() = delete;

    /**
     * @param vname Name of the virus.
     * @param n Number of agents.
     * @param prevalence Initial share of infected (exposed) agents.
     * @param transmission_prob Probability of transmission per contact.
     * @param latent_period Mean latent period (days; geometric).
     * @param recovery_rate Daily probability of recovery while infected.
     * @param hospitalization_rate Daily probability of becoming severe.
     * @param hospitalization_period Mean days severe (in hospital).
     * @param detection_delay Mean days until an infected agent is detected
     * (negative: never; 0: at once).
     * @param quarantine_period Days of quarantine (negative: no quarantine).
     * @param quarantine_willingness Share of agents willing to quarantine.
     * @param isolation_willingness Share of agents willing to isolate.
     * @param isolation_period Days of isolation from the onset (negative: no
     * isolation).
     * @param contact_tracing_prob Probability that a contact is traced.
     * @param contact_tracing_window Days back the contacts are traced.
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

    void reset() override;
    std::unique_ptr< Model > clone_ptr() override;

    std::vector< size_t > get_agent_quarantine_triggered() const
    {
        return agent_quarantine_triggered;
    };

    std::vector< bool > get_quarantine_willingness() const
    {
        return quarantine_willingness;
    };

    std::vector< bool > get_isolation_willingness() const
    {
        return isolation_willingness;
    };


    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};

inline void ModelSEIRQuarantine::reset()
{

    // The quarantine process reads the contacts of every layer
    for (size_t i = 0u; i < this->get_n_layers(); ++i)
        if (!this->get_layer(i).is_tracing())
            this->get_layer(i).enable_tracing();

    Model::reset();

    // Who is willing to quarantine, and to isolate (drawn once per run)
    quarantine_willingness.assign(this->size(), false);
    isolation_willingness.assign(this->size(), false);
    for (size_t idx = 0; idx < quarantine_willingness.size(); ++idx)
    {
        quarantine_willingness[idx] =
            this->runif() < EPI_PAR(this, "Quarantine willingness");
        isolation_willingness[idx] =
            this->runif() < EPI_PAR(this, "Isolation willingness");
    }

    agent_quarantine_triggered.assign(this->size(), 0u);

    return;

}

inline std::unique_ptr<Model> ModelSEIRQuarantine::clone_ptr()
{
    return std::make_unique<ModelSEIRQuarantine>(*this);
}

inline void ModelSEIRQuarantine::_quarantine_process(Model * m)
{

    auto * model = model_cast<ModelSEIRQuarantine>(m);
    const int isolated = m->transition_id("Isolated");
    const int quarantined = m->transition_id("Quarantined");
    const bool isolation_on = EPI_PAR(m, "Isolation period") >= 0;

    // 1. Detection of the free infected agents
    const epiworld_double delay = EPI_PAR(m, "Detection delay");
    if (delay >= 0.0)
    {

        for (auto & p : m->get_agents())
        {

            // From the day after the onset: its contacts as an infectious
            // agent are in the logs by then
            const size_t id = static_cast< size_t >(p.get_id());
            if (
                (m->get_stage(p) != INFECTED) || (m->get_status(p) != FREE) ||
                (model->agent_quarantine_triggered[id] != QUARANTINE_PROCESS_INACTIVE) ||
                (m->get_stage_entered(id, INFECTED) >= m->today())
            )
                continue;

            if ((delay > 0.0) && (m->runif() >= 1.0 / delay))
                continue;

            model->agent_quarantine_triggered[p.get_id()] = QUARANTINE_PROCESS_ACTIVE;
            if (isolation_on && model->isolation_willingness[p.get_id()])
                m->fire(isolated, p);

        }

    }

    // 2. Contact tracing from the detected agents
    for (size_t agent_i = 0u; agent_i < m->size(); ++agent_i)
    {

        if (model->agent_quarantine_triggered[agent_i] != QUARANTINE_PROCESS_ACTIVE)
            continue;

        model->agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;

        if (EPI_PAR(m, "Quarantine period") < 0)
            continue;

        const auto success_rate = EPI_PAR(m, "Contact tracing probability");
        const auto days_prior = EPI_PAR(m, "Contact tracing window");

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

                if (m->runif() > success_rate)
                    continue;

                auto [contact_id, contact_date] = ct.get_contact(agent_i, contact_i);

                if ((static_cast< double >(m->today()) -
                    static_cast< double >(contact_date)) > days_prior)
                    continue;

                auto & agent = m->get_agent(contact_id);
                const int s = m->get_stage(agent);
                if ((m->get_status(agent) != FREE) || (s > INFECTED))
                    continue;

                // Agents with some tool are not quarantined
                if (agent.get_n_tools() != 0u)
                    continue;

                if (!model->quarantine_willingness[contact_id])
                    continue;

                if (s != INFECTED)
                    m->fire(quarantined, agent);
                else if (isolation_on && model->isolation_willingness[contact_id])
                    m->fire(isolated, agent);

            }

        }

    }

}

inline ModelSEIRQuarantine::ModelSEIRQuarantine(
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

    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(latent_period, "Latent period");
    this->add_param(hospitalization_rate, "Hospitalization rate");
    this->add_param(hospitalization_period, "Hospitalization period");
    this->add_param(detection_delay, "Detection delay");
    this->add_param(quarantine_period, "Quarantine period");
    this->add_param(quarantine_willingness, "Quarantine willingness");
    this->add_param(isolation_willingness, "Isolation willingness");
    this->add_param(isolation_period, "Isolation period");
    this->add_param(contact_tracing_prob, "Contact tracing probability");
    this->add_param(contact_tracing_window, "Contact tracing window");

    // The disease
    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Exposed", stage::carrier);
    virus.add_stage("Infected", stage::infectious);
    virus.add_stage("Severe", stage::carrier);
    virus.add_stage("Recovered");
    virus.infection(SUSCEPTIBLE, EXPOSED, "Transmission probability");
    virus.after(EXPOSED, INFECTED, Duration::geometric("Latent period"));
    virus.rate(INFECTED, RECOVERED, "Recovery rate").effect(Effect::recovery);
    virus.rate(INFECTED, SEVERE, "Hospitalization rate");
    const int out_severe =
        virus.rate(SEVERE, RECOVERED, Value::inverse("Hospitalization period"));
    this->add_virus(virus, assign::share(prevalence));

    // Where agents are
    Statuses s;
    s.add_status("Free");
    s.add_status("Quarantined");
    s.add_status("Isolated");
    s.add_status("Hospitalized");
    s.set_contact(QUARANTINED, 0.0);
    s.set_contact(ISOLATED, 0.0);
    s.set_contact(HOSPITALIZED, 0.0);

    s.external(FREE, QUARANTINED, "Quarantined");
    s.after(QUARANTINED, FREE, Duration::fixed("Quarantine period"))
        .name("Quarantine ends");
    s.external(FREE, ISOLATED, "Isolated");
    const int onset_q = s.external({FREE, QUARANTINED}, ISOLATED, "Onset in quarantine");
    s.after(
        ISOLATED, FREE,
        Duration::fixed("Isolation period").since(INFECTED)
    ).name("Isolation ends");
    const int admit = s.external({FREE, QUARANTINED, ISOLATED}, HOSPITALIZED, "Admitted");
    const int discharge = s.external(HOSPITALIZED, FREE, "Discharged");
    this->set_statuses(std::move(s));

    // How the two are linked
    this->react(Trigger::entered(INFECTED).was({QUARANTINED}), onset_q);
    this->react(Trigger::entered(SEVERE), admit);
    this->react(
        Trigger::entered(SEVERE),
        [](Agent & p, Model & m) { m.record_hospitalization(p); }
    );
    this->react(Trigger::fired(out_severe), discharge);

    this->add_globalevent(_quarantine_process, "Quarantine process");

    this->set_name("SEIR with Quarantine");

}

inline ModelInfo ModelSEIRQuarantine::model_info()
{

    using uint_t = epiworld_fast_uint;
    using int_t  = epiworld_fast_int;
    const double inf = std::numeric_limits< double >::infinity();

    return {"SEIRQuarantine",
        "SEIR with hospitalization, isolation and quarantine",
        "SEIR where detected cases are isolated and their contacts, traced "
        "through the layers that trace, quarantined.",
        "Disease",
        {
            ParamInfo::agents(10000),
            ParamInfo::prevalence(0.01),
            ParamInfo::probability("transmission_prob", "Transmission probability",
                0.05, "Probability that a contact with an infected agent infects."),
            ParamInfo::duration("latent_period", "Latent period", 7.0,
                "Mean days in the exposed state (infected, not yet infectious)."),
            ParamInfo::probability("recovery_rate", "Recovery rate", 0.14,
                "Daily probability that an infected agent recovers."),
            ParamInfo::probability("hospitalization_rate", "Hospitalization rate",
                0.05, "Daily probability that an infected agent is hospitalized."),
            ParamInfo::duration("hospitalization_period", "Hospitalization period",
                7.0, "Average number of days in the hospital."),
            {"detection_delay", "Detection delay", 3.0, -1.0, inf, false,
                "Mean days before an infected agent is detected (0: at once; "
                "negative: never)."},
            ParamInfo::period("quarantine_period", "Quarantine period", 14.0,
                "Days in quarantine for the contacts of a detected case "
                "(negative: no quarantine)."),
            ParamInfo::probability("quarantine_willingness",
                "Quarantine willingness", 0.8,
                "Probability that a contact complies with quarantine."),
            ParamInfo::probability("isolation_willingness",
                "Isolation willingness", 0.5,
                "Probability that a detected case complies with isolation."),
            ParamInfo::period("isolation_period", "Isolation period", 7.0,
                "Days in isolation for a detected case (negative: no "
                "isolation)."),
            ParamInfo::probability("contact_tracing_prob",
                "Contact tracing probability", 1.0,
                "Probability that a contact is traced."),
            ParamInfo::count("contact_tracing_window", "Contact tracing window",
                4.0, "Days before detection whose contacts are traced.")
        },
        [](const ModelArgs & a) -> std::unique_ptr< Model > {
            return std::make_unique< ModelSEIRQuarantine >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("latent_period"), a("recovery_rate"),
                a("hospitalization_rate"), a("hospitalization_period"),
                a("detection_delay"),
                static_cast< int_t >(a("quarantine_period")),
                a("quarantine_willingness"), a("isolation_willingness"),
                static_cast< int_t >(a("isolation_period")),
                a("contact_tracing_prob"),
                static_cast< uint_t >(a("contact_tracing_window"))
            );
        }};
}

#endif
