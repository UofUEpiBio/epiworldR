
#ifndef MEASLESSCHOOL_HPP
#define MEASLESSCHOOL_HPP

#include <cassert>
#include "../tools/vaccine.hpp"
#include "../model-bones.hpp"
#include "../globalevents/pep.hpp"

#define LOCAL_UPDATE_FUN(name) \
    template<typename TSeq> \
    inline void ModelMeaslesSchool<TSeq>:: name \
    (Agent<TSeq> * p, Model<TSeq> * m)

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
 * 
 * ![Model Diagram](../assets/img/measlesschool.png)
 * 
 * In the case of post-exposure prophylaxis (PEP), agents who are quarantined
 * and exposed or susceptible can receive PEP with a certain willingness and
 * efficacy. If they receive PEP, then they have a probability of moving to the
 * recovered state, which is a function of the PEP efficacy.
 * 
 * 
 * @ingroup disease_specific
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelMeaslesSchool final : public Model<TSeq> {

private:

    /**
     * @brief The function that updates the model.
     */
    ///@{
    static void _update_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_rash(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_isolated(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_isolated_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_q_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_q_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_q_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_q_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void _update_hospitalized(Agent<TSeq> * p, Model<TSeq> * m);
    ///@}

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
    static void _quarantine_agents(Model<TSeq> * m);
    
    // Update which agents are infectious for contact
    void _update_infectious();

public:

    /**
     * @name Model States
     * @brief The different states of the model.
     */
    ///@{
    static constexpr epiworld_fast_uint SUSCEPTIBLE             = 0u;
    static constexpr epiworld_fast_uint EXPOSED                 = 1u;
    static constexpr epiworld_fast_uint PRODROMAL               = 2u;
    static constexpr epiworld_fast_uint RASH                    = 3u;
    static constexpr epiworld_fast_uint ISOLATED                = 4u;
    static constexpr epiworld_fast_uint ISOLATED_RECOVERED      = 5u;
    static constexpr epiworld_fast_uint DETECTED_HOSPITALIZED   = 6u;
    static constexpr epiworld_fast_uint QUARANTINED_EXPOSED     = 7u;
    static constexpr epiworld_fast_uint QUARANTINED_SUSCEPTIBLE = 8u;
    static constexpr epiworld_fast_uint QUARANTINED_PRODROMAL   = 9u;
    static constexpr epiworld_fast_uint QUARANTINED_RECOVERED   = 10u;
    static constexpr epiworld_fast_uint HOSPITALIZED            = 11u;
    static constexpr epiworld_fast_uint RECOVERED               = 12u;
    ///@}
    
    // Default constructor
    ModelMeaslesSchool() {};

    ModelMeaslesSchool(
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
        epiworld_fast_int isolation_period,
        // Policy parameters 2
        epiworld_double pep_efficacy = 0.0,
        epiworld_double pep_willingness = 0.0
    );
    ///@}

    std::vector<Agent<TSeq> *> infectious; ///< Agents infectious for contact

    bool system_quarantine_triggered = false;

    std::vector< int > day_flagged; ///< Either detected or started quarantine
    std::vector< int > day_rash_onset; ///< Day of rash onset
    std::vector< int > has_pep;

    void reset() override;

    std::unique_ptr< Model<TSeq> > clone_ptr() override;
    void next() override;

};

template<typename TSeq>
inline void ModelMeaslesSchool<TSeq>::_quarantine_agents(Model<TSeq> * m) {

    auto * model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);

    // Iterating through the new cases
    if (!model->system_quarantine_triggered)
        return;

    // Quarantine and isolation can be shut off if negative
    if (
        (model->par("Quarantine period") < 0) &&
        (model->par("Isolation period") < 0)
    )
        return;

    // Capturing the days that matter and the probability of success
    epiworld_double willingness = model->par("Quarantine willingness");

    // Iterating through the
    for (size_t i = 0u; i < model->size(); ++i) {

        auto & agent = model->get_agent(i);
        auto agent_state = agent.get_state();

        // Already quarantined or isolated
        if (agent_state >= RASH)
            continue;

        // If the agent has a vaccine, then no need for quarantine
        if (agent.get_n_tools() != 0u)
            continue;


        // Quarantine will depend on the willingness of the agent
        // to be quarantined. If negative, then quarantine never happens.
        if (
            (model->par("Quarantine period") >= 0) &&
            (model->runif() < willingness)
        )
        {

            if (agent_state == SUSCEPTIBLE)
                agent.change_state(*model, QUARANTINED_SUSCEPTIBLE);
            else if (agent_state == EXPOSED)
                agent.change_state(*model, QUARANTINED_EXPOSED);
            else if (agent_state == PRODROMAL)
                agent.change_state(*model, QUARANTINED_PRODROMAL);

            // And we add the day of quarantine
            model->day_flagged[i] = model->today();

        }

    }

    // Setting the quarantine process off
    model->system_quarantine_triggered = false;

    return;

}

template<typename TSeq>
inline void ModelMeaslesSchool<TSeq>::reset() {

    Model<TSeq>::reset();

    this->system_quarantine_triggered = false;

    this->day_flagged.assign(this->size(), 0);
    this->day_rash_onset.assign(this->size(), 0);
    this->has_pep.assign(this->size(), false);

    this->_update_infectious();
    return;

}

template<typename TSeq>
inline void ModelMeaslesSchool<TSeq>::_update_infectious() {

    #ifdef EPI_DEBUG
    // All agents with state >= EXPOSED should have a virus
    for (auto & agent: this->get_agents())
    {
        auto s = agent.get_state();
        static const std::vector< int > states_with_virus = {
            EXPOSED, PRODROMAL, RASH, ISOLATED, DETECTED_HOSPITALIZED,
            QUARANTINED_EXPOSED, QUARANTINED_PRODROMAL, HOSPITALIZED
        };

        if (IN(s, states_with_virus))
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

        if ((s < RASH) || (s == RECOVERED))
            ++n_available;

    }

    // Assumes fixed contact rate throughout the simulation
    // but corrects for the number of available agents.
    double p_contact = 0.0;
    if (n_available > 0)
    {
        p_contact = this->par("Contact rate")/
            static_cast< epiworld_double >(n_available);
    }

    // Notice this is for sampling with replacement
    // from the list of infected individuals.
    // This is a partial drawing which, complemented with
    // drawing from the non-infected individuals, yields
    // a Binomial(n, contact_rate/n) number of contacts.
    this->set_rand_binom(
        static_cast<int>(this->infectious.size()),
        p_contact > 1.0 ? 1.0 : p_contact
    );

}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelMeaslesSchool<TSeq>::clone_ptr()
{

    return std::make_unique<ModelMeaslesSchool<TSeq>>(*this);

}

LOCAL_UPDATE_FUN(_update_susceptible) {

    // How many contacts to draw
    int ndraw = m->rbinom();

    if (ndraw == 0)
        return;

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);
    size_t n_infectious = model->infectious.size();

    if (n_infectious == 0)
        return;

    // Drawing from the set
    int nviruses_tmp = 0;
    int i = 0;
    auto & _ref = *m;
    while (i < ndraw)
    {
        // Picking the actual contacts
        int which = m->runif_int(0, n_infectious - 1);

        Agent<> & neighbor = *model->infectious[which];

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
            (1.0 - p->get_susceptibility_reduction(v, _ref)) *
            v->get_prob_infecting(m) *
            (1.0 - neighbor.get_transmission_reduction(v, _ref))
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

    p->set_virus(*m, *m->array_virus_tmp[which]);

    return;

};

LOCAL_UPDATE_FUN(_update_exposed) {

    // if (InterventionPEP<TSeq>::agent_recovers(*p, *m, RECOVERED))
    //     return;

    if (m->runif() < (1.0/p->get_virus()->get_incubation(m)))
        p->change_state(*m, ModelMeaslesSchool<TSeq>::PRODROMAL);

    return;

};

LOCAL_UPDATE_FUN(_update_prodromal) {

    if (m->runif() < (1.0/m->par("Prodromal period")))
    {

        auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);
        model->day_rash_onset[p->get_id()] = m->today();
        p->change_state(*m, ModelMeaslesSchool<TSeq>::RASH);

    }

    return;

};

LOCAL_UPDATE_FUN(_update_rash) {


    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);

    #ifdef EPI_DEBUG
    if (static_cast<int>(model->day_flagged.size()) <= p->get_id())
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
    auto which = m->sample_from_probs(2);

    // Recovers (which == 0 fires with probability 1/rash_period)
    if (which == 0)
    {
        p->rm_virus(*m, detected ? ISOLATED_RECOVERED: RECOVERED);
    }
    else if (which == 1)
    {
        // If hospitalized, then the agent is removed from the system
        // effectively
        model->record_hospitalization(*p);
        p->change_state(*m, detected ? DETECTED_HOSPITALIZED : HOSPITALIZED);
    }
    else if (which > 2)
    {
        throw std::logic_error("The roulette returned an unexpected value.");
    }
    else if (detected)
    {
        // If the agent is not hospitalized, then it is moved to
        // isolation.
        p->change_state(*m, ISOLATED);
    }

};

LOCAL_UPDATE_FUN(_update_isolated) {

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);

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
    auto which = m->sample_from_probs(2);

    // Recovers (which == 0 fires with probability 1/rash_period)
    if (which == 0u)
    {
        p->rm_virus(*m, unisolate? RECOVERED: ISOLATED_RECOVERED);
    }

    // If hospitalized, then the agent is removed from the system
    else if (which == 1u)
    {
        model->record_hospitalization(*p);
        p->change_state(*m, unisolate? HOSPITALIZED: DETECTED_HOSPITALIZED);
    }
    // If neither hospitalized nor recovered, then the agent is
    // still under isolation, unless the quarantine period is over.
    else if (unisolate)
    {
        p->change_state(*m, RASH);
    }

}

LOCAL_UPDATE_FUN(_update_isolated_recovered) {

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_rash_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    if (unisolate)
        p->change_state(*m, RECOVERED);

}

LOCAL_UPDATE_FUN(_update_q_exposed) {

    #ifdef EPI_DEBUG
    if (m->par("PEP willingness") >= 0.999999999)
    {
        throw std::logic_error(
            std::string("This shouldn't happen. ") +
            std::string("When PEP willingness is 1, then no agent should be ") +
            std::string("in QUARANTINED_EXPOSED state.")
        );
    }
    #endif

    // How many days since quarantine started
    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);
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
        p->change_state(*m, unquarantine ? PRODROMAL : QUARANTINED_PRODROMAL);

    }
    else if (unquarantine)
    {
        p->change_state(*m, 
            EXPOSED
        );
    }

}

LOCAL_UPDATE_FUN(_update_q_susceptible) {

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);
    int days_since =
        m->today() - model->day_flagged[p->get_id()];

    if (days_since >= m->par("Quarantine period"))
        p->change_state(*m, SUSCEPTIBLE);

}

LOCAL_UPDATE_FUN(_update_q_prodromal) {

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);

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
        p->change_state(*m, ISOLATED);
    }
    else
    {

        if (unquarantine)
            p->change_state(*m, PRODROMAL);

    }

}

LOCAL_UPDATE_FUN(_update_q_recovered) {

    auto* model = model_cast<ModelMeaslesSchool<TSeq>,TSeq>(m);
    int days_since = m->today() - model->day_flagged[p->get_id()];

    if (days_since >= m->par("Quarantine period"))
        p->change_state(*m, RECOVERED);

}

LOCAL_UPDATE_FUN(_update_hospitalized) {

    // The agent is removed from the system
    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_virus(*m, RECOVERED);

    return;

}


template<typename TSeq>
inline ModelMeaslesSchool<TSeq>::ModelMeaslesSchool(
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
    epiworld_fast_int isolation_period,
    epiworld_double pep_efficacy,
    epiworld_double pep_willingness
) {

    this->add_state("Susceptible",             this->_update_susceptible);
    this->add_state("Exposed",                 this->_update_exposed);
    this->add_state("Prodromal",               this->_update_prodromal);
    this->add_state("Rash",                    this->_update_rash);
    this->add_state("Isolated",                this->_update_isolated);
    this->add_state("Isolated Recovered",      this->_update_isolated_recovered);
    this->add_state("Detected Hospitalized",   this->_update_hospitalized);
    this->add_state("Quarantined Exposed",     this->_update_q_exposed);
    this->add_state("Quarantined Susceptible", this->_update_q_susceptible);
    this->add_state("Quarantined Prodromal",   this->_update_q_prodromal);
    this->add_state("Quarantined Recovered",   this->_update_q_recovered);
    this->add_state("Hospitalized",            this->_update_hospitalized);
    this->add_state("Recovered");

    // Adding the model parameters
    this->add_param(contact_rate, "Contact rate");
    this->add_param(transmission_rate, "Transmission rate");
    this->add_param(incubation_period, "Incubation period");
    this->add_param(prodromal_period, "Prodromal period");
    this->add_param(rash_period, "Rash period");
    this->add_param(days_undetected, "Days undetected");
    this->add_param(quarantine_period, "Quarantine period");
    this->add_param(quarantine_willingness, "Quarantine willingness");
    this->add_param(isolation_period, "Isolation period");
    this->add_param(hospitalization_rate, "Hospitalization rate");
    this->add_param(hospitalization_period, "Hospitalization period");
    this->add_param(prop_vaccinated, "Vaccination rate");
    this->add_param(vax_efficacy, "Vax efficacy");
    this->add_param(vax_reduction_recovery_rate, "(IGNORED) Vax improved recovery");
    this->add_param(pep_efficacy, "PEP efficacy");
    this->add_param(pep_willingness, "PEP willingness");

    // Designing the disease
    Virus<> measles("Measles");
    measles.set_state(EXPOSED, RECOVERED);
    measles.set_prob_infecting("Transmission rate");
    measles.set_prob_recovery("Rash period");
    measles.set_incubation("Incubation period");
    measles.set_distribution(
        distribute_virus_randomly(n_exposed, false)
    );

    this->add_virus(measles);

    // Designing the vaccine
    ToolVaccine<TSeq> vaccine("MMR");
    vaccine.set_susceptibility_reduction(this->par("Vax efficacy"));
    vaccine.set_distribution(distribute_tool_randomly(prop_vaccinated, true));
    this->add_tool(vaccine);

    this->queuing_off();

    // Quarantine process will be automatically triggered
    // at the end of the day
    auto quarantine_event = GlobalEvent<TSeq>(
        this->_quarantine_agents, "Quarantine process"
    );
    this->add_globalevent(quarantine_event);

    // Creating the PEP intervention and 
    // setting it up so we can call it as a global event.
    InterventionPEP<TSeq> pep(
        "PEP efficacy",
        "PEP willingness",
        {QUARANTINED_EXPOSED, QUARANTINED_SUSCEPTIBLE},
        {RECOVERED, SUSCEPTIBLE}
    );

    pep.set_name("PEP intervention");

    this->add_globalevent(pep);

    // Setting the population
    this->agents_empty_graph(n);

}

template<typename TSeq>
inline void ModelMeaslesSchool<TSeq>::next() {

    this->_update_infectious();
    Model<TSeq>::next();

}

#undef LOCAL_UPDATE_FUN
#endif
