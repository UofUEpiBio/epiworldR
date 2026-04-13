#ifndef EPIWORLD_MODELS_MEASLESMIXING_HPP
#define EPIWORLD_MODELS_MEASLESMIXING_HPP

#include "../tools/vaccine.hpp"
#include "../model-bones.hpp"

#define MM(i, j, n) \
    j * n + i

/**
 * @file measlesmixing.hpp
 * @brief Template for a Measles model with population mixing, quarantine, and contact tracing
 */

/**
 * @brief Measles model with population mixing, quarantine, and contact tracing
 *
 * This class implements a Measles epidemiological model based on the SEIR framework
 * with additional features including:
 * - Population mixing based on contact matrices
 * - Measles-specific disease progression: Susceptible → Latent → Prodromal → Rash
 * - Prodromal individuals are infectious (replace the "Infected" state in SEIR)
 * - Rash individuals are no longer infectious but can be detected for isolation
 * - Quarantine measures for latent contacts during contact tracing
 * - Isolation policies for detected individuals during the rash state
 * - Contact tracing with configurable success rates
 * - Hospitalization of severe cases
 * - Individual willingness to comply with public health measures
 *
 * The model supports 13 distinct states:
 * - Susceptible: Individuals who can become infected
 * - Latent: Infected but not yet infectious (incubation period)
 * - Prodromal: Infectious individuals in the community (replaces "Infected" in SEIR)
 * - Rash: Non-infectious individuals with visible symptoms (detection occurs here)
 * - Isolated: Detected individuals in self-isolation
 * - Isolated Recovered: Recovered individuals still in isolation
 * - Detected Hospitalized: Hospitalized individuals who were contact-traced
 * - Quarantined Latent: Latent individuals in quarantine due to contact tracing
 * - Quarantined Susceptible: Susceptible individuals in quarantine due to contact tracing
 * - Quarantined Prodromal: Prodromal individuals in quarantine due to contact tracing
 * - Quarantined Recovered: Recovered individuals in quarantine due to contact tracing
 * - Hospitalized: Individuals requiring hospital care
 * - Recovered: Individuals who have recovered and gained immunity
 *
 * ![Model Diagram](../assets/img/measlesmixing.png)
 *
 * @tparam TSeq Type for genetic sequences (default: EPI_DEFAULT_TSEQ)
 * @ingroup disease_specific
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelMeaslesMixing final : public Model<TSeq>
{
private:

    // Vector of vectors of infected agents (prodromal agents are infectious)
    std::vector< size_t > infectious;

    // Number of infectious agents in each group
    std::vector< size_t > n_infectious_per_group;

    // Where the agents start in the `infectious` vector
    std::vector< size_t > entity_indices;

    void m_update_infectious_list();
    std::vector< size_t > sampled_agents;
    size_t sample_agents(
        Agent<TSeq> * agent,
        std::vector< size_t > & sampled_agents
        );
    std::vector< double > adjusted_contact_rate;
    std::vector< double > contact_matrix;

    #ifdef EPI_DEBUG
    std::vector< int > sampled_sizes;
    #endif

    // Update functions
    static void m_update_susceptible(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_latent(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_rash(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_suscep(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_latent(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_prodromal(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_recovered(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_hospitalized(Agent<TSeq> * p, Model<TSeq> * m);

    // Data about the quarantine process
    std::vector< bool > quarantine_willingness; ///< Indicator for quarantine willingness
    std::vector< bool > isolation_willingness; ///< Indicator for isolation willingness
    std::vector< size_t > agent_quarantine_triggered; ///< Whether the quarantine process has started
    std::vector< int > day_flagged; ///< Either detected or started quarantine
    std::vector< int > day_rash_onset; ///< Day of rash onset
    std::vector< int > day_latent; ///< Day of latent infection

    void m_quarantine_process();
    void m_update_model();

    // We will limit tracking to up to EPI_MAX_TRACKING
    ContactTracing contact_tracing;

public:

    static const int SUSCEPTIBLE              = 0;
    static const int LATENT                   = 1;
    static const int PRODROMAL                = 2;
    static const int RASH                     = 3;
    static const int ISOLATED                 = 4;
    static const int ISOLATED_RECOVERED       = 5;
    static const int DETECTED_HOSPITALIZED    = 6;
    static const int QUARANTINED_LATENT       = 7;
    static const int QUARANTINED_SUSCEPTIBLE  = 8;
    static const int QUARANTINED_PRODROMAL    = 9;
    static const int QUARANTINED_RECOVERED    = 10;
    static const int HOSPITALIZED             = 11;
    static const int RECOVERED                = 12;

    static const size_t QUARANTINE_PROCESS_INACTIVE = 0u;
    static const size_t QUARANTINE_PROCESS_ACTIVE   = 1u;
    static const size_t QUARANTINE_PROCESS_DONE     = 2u;

    ModelMeaslesMixing() {};

    /**
     * @brief Constructs a ModelMeaslesMixing object.
     *
     * @param n The number of entities in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param contact_rate The contact rate between entities in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param vax_efficacy The efficacy of the vaccine.
     * @param vax_reduction_recovery_rate The reduction in recovery rate due to the vaccine.
     * @param incubation_period The incubation period of the disease in the model.
     * @param prodromal_period The prodromal period of the disease in the model.
     * @param rash_period The rash period of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model. Specified in
     * column-major order.
     * @param hospitalization_rate The rate at which infected individuals are hospitalized.
     * @param hospitalization_period The average duration of hospitalization in days.
     * @param days_undetected The average number of days an infected individual remains undetected.
     * @param quarantine_period The duration of quarantine in days for latent contacts.
     * @param quarantine_willingness The proportion of individuals willing to comply with quarantine measures.
     * @param isolation_willingness The proportion of individuals willing to self-isolate when detected.
     * @param isolation_period The duration of isolation in days for detected infected individuals.
     * @param prop_vaccinated The proportion of vaccinated agents.
     * @param contact_tracing_success_rate The probability of successfully identifying and tracing contacts (default: 1.0).
     * @param contact_tracing_days_prior The number of days prior to detection for which contacts are traced (default: 4).
     */
    ModelMeaslesMixing(
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double vax_efficacy,
        epiworld_double vax_reduction_recovery_rate,
        epiworld_double incubation_period,
        epiworld_double prodromal_period,
        epiworld_double rash_period,
        std::vector< double > contact_matrix,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double days_undetected,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_double isolation_willingness,
        epiworld_fast_int isolation_period,
        epiworld_double prop_vaccinated,
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
     * - [0]: The proportion of initially infected individuals who start in the latent state.
     * - [1]: The proportion of initially non-infected individuals who have recovered (immune).
     * @param queue_ Optional vector for queuing specifications (default: empty).
     */
    ModelMeaslesMixing<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    ) override;

    /**
     * @brief Set the contact matrix for population mixing
     * @param cmat Contact matrix specifying interaction rates between groups
     */
    void set_contact_matrix(std::vector< double > cmat)
    {
        contact_matrix = cmat;
        return;
    };

    /**
     * @brief Get the current contact matrix
     * @return Vector representing the contact matrix
     */
    std::vector< double > get_contact_matrix() const
    {
        return contact_matrix;
    };

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

    // Overriding the next() function to include the model update
    void next() override;

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_model()
{
    this->m_quarantine_process();
    this->events_run();
    this->m_update_infectious_list();

}

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_infectious_list()
{

    auto & agents = this->get_agents();

    // Resetting infectious list
    std::fill(n_infectious_per_group.begin(), n_infectious_per_group.end(), 0u);

    // Resetting the number of available contacts
    adjusted_contact_rate.assign(this->entities.size(), 0.0);

    for (const auto & a : agents)
    {

        if (a.get_state() == PRODROMAL)
        {
            if (a.get_n_entities() > 0u)
            {
                const auto & entity = a.get_entity(0u, *this);
                infectious[
                    // Position of the group in the `infectious` vector
                    entity_indices[entity.get_id()] +
                    // Position of the agent in the group
                    n_infectious_per_group[entity.get_id()]++
                ] = a.get_id();

            }
        }

        // Setting how many agents are available for contact
        if (
            ((a.get_state() < RASH) || (a.get_state() == RECOVERED)) &&
            (a.get_n_entities() > 0u)
        )
        {
            adjusted_contact_rate[
                a.get_entity(0u, *this).get_id()
            ] += 1.0;
        }

    }

    // This simplifies calculations later
    for (auto & rate: adjusted_contact_rate)
    {
        if (rate > 0.0)
            rate = this->get_param("Contact rate") / rate;
        else
            rate = 0.0;  // No available contacts in this group

        if (rate > 1.0)
            rate = 1.0;
    }

    return;

}

template<typename TSeq>
inline size_t ModelMeaslesMixing<TSeq>::sample_agents(
    Agent<TSeq> * agent,
    std::vector< size_t > & sampled_agents
    )
{

    size_t agent_group_id = agent->get_entity(0u, *this).get_id();
    size_t ngroups = this->entities.size();

    int samp_id = 0;
    for (size_t g = 0; g < ngroups; ++g)
    {

        size_t group_size = n_infectious_per_group[g];

        if (group_size == 0u)
            continue;

        // How many from this entity?
        int nsamples = this->rbinom(
            group_size,
            adjusted_contact_rate[g] * contact_matrix[
                MM(agent_group_id, g, ngroups)
            ]
        );

        if (nsamples == 0)
            continue;

        // Sampling from the entity
        for (int s = 0; s < nsamples; ++s)
        {

            // Randomly selecting an agent
            int which = this->runif() * group_size;

            // Correcting overflow error
            if (which >= static_cast<int>(group_size))
                which = static_cast<int>(group_size) - 1;

            #ifdef EPI_DEBUG
            auto & a = this->population.at(infectious.at(entity_indices[g] + which));
            #else
            auto & a = this->get_agent(infectious[entity_indices[g] + which]);
            #endif

            #ifdef EPI_DEBUG
            if (a.get_state() != PRODROMAL)
                throw std::logic_error(
                    "The agent is not in prodromal state, but it should be."
                );
            #endif

            // Can't sample itself
            if (a.get_id() == agent->get_id())
                continue;

            sampled_agents[samp_id++] = a.get_id();

        }

    }

    return samp_id;

}

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::reset()
{

    Model<TSeq>::reset();

    // Checking contact matrix's rows add to one
    size_t nentities = this->entities.size();
    if (this->contact_matrix.size() !=  nentities*nentities)
        throw std::length_error(
            std::string("The contact matrix must be a square matrix of size ") +
            std::string("nentities x nentities. ") +
            std::to_string(this->contact_matrix.size()) +
            std::string(" != ") + std::to_string(nentities*nentities) +
            std::string(".")
            );

    for (size_t i = 0u; i < this->entities.size(); ++i)
    {
        double sum = 0.0;
        for (size_t j = 0u; j < this->entities.size(); ++j)
        {
            if (this->contact_matrix[MM(i, j, nentities)] < 0.0)
                throw std::range_error(
                    std::string("The contact matrix must be non-negative. ") +
                    std::to_string(this->contact_matrix[MM(i, j, nentities)]) +
                    std::string(" < 0.")
                    );
            sum += this->contact_matrix[MM(i, j, nentities)];
        }
        if (sum < 0.999 || sum > 1.001)
            throw std::range_error(
                std::string("The contact matrix must have rows that add to one. ") +
                std::to_string(sum) +
                std::string(" != 1.")
                );
    }

    // Do it the first time only
    sampled_agents.resize(this->size());

    // We only do it once
    n_infectious_per_group.assign(this->entities.size(), 0u);

    // We are assuming one agent per entity
    infectious.assign(this->size(), 0u);

    // This will say when do the groups start in the `infectious` vector
    entity_indices.assign(this->entities.size(), 0u);

    for (size_t i = 1u; i < this->entities.size(); ++i)
    {

        entity_indices[i] +=
            this->entities[i - 1].size() +
            entity_indices[i - 1]
            ;

    }

    this->m_update_infectious_list();

    // Setting up the quarantine parameters
    quarantine_willingness.assign(this->size(), false);
    isolation_willingness.assign(this->size(), false);
    for (size_t idx = 0; idx < quarantine_willingness.size(); ++idx)
    {
        quarantine_willingness[idx] =
            this->runif() < this->par("Quarantine willingness");
        isolation_willingness[idx] =
            this->runif() < this->par("Isolation willingness");
    }

    agent_quarantine_triggered.assign(this->size(), 0u);
    day_flagged.assign(this->size(), 0);
    day_rash_onset.assign(this->size(), 0);
    day_latent.assign(this->size(), 0);

    // Contact tracing
    contact_tracing.reset(this->size(), EPI_MAX_TRACKING);

    return;

}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelMeaslesMixing<TSeq>::clone_ptr()
{

    return std::make_unique<ModelMeaslesMixing<TSeq>>(*this);

}

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_susceptible(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    if (p->get_n_entities() == 0)
        return;

    // Downcasting to retrieve the sampler attached to the
    // class
    auto * m_down = model_cast<ModelMeaslesMixing<TSeq>, TSeq>(m);

    size_t ndraws = m_down->sample_agents(p, m_down->sampled_agents);

    #ifdef EPI_DEBUG
    m_down->sampled_sizes.push_back(static_cast<int>(ndraws));
    #endif

    if (ndraws == 0u)
        return;

    // Drawing from the set
    int nviruses_tmp = 0;
    auto & m_ref = *m;
    for (size_t n = 0u; n < ndraws; ++n)
    {

        auto & neighbor = m->get_agent(m_down->sampled_agents[n]);

        auto & v = neighbor.get_virus();

        #ifdef EPI_DEBUG
        if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
            throw std::logic_error(
                "Trying to add an extra element to a temporal array outside of the range."
            );
        #endif

        // Adding the current agent to the tracked interactions
        m_down->contact_tracing.add_contact(neighbor.get_id(), p->get_id(), m->today());

        /* And it is a function of susceptibility_reduction as well */
        m->array_double_tmp[nviruses_tmp] =
            (1.0 - p->get_susceptibility_reduction(v, m_ref)) *
            v->get_prob_infecting(m) *
            (1.0 - neighbor.get_transmission_reduction(v, m_ref))
            ;

        m->array_virus_tmp[nviruses_tmp++] = &(*v);

    }

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return;

    p->set_virus(*m, *m->array_virus_tmp[which], LATENT);

    return;

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_latent(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // Getting the virus
    auto & v = p->get_virus();

    // Does the agent become prodromal (infectious)?
    if (m->runif() < 1.0/(v->get_incubation(m)))
    {

        p->change_state(*m, PRODROMAL);

        return;

    }

    return;

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_prodromal(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>, TSeq>(m);

    // Does the agent transition to rash?
    if (m->runif() < 1.0/m->par("Prodromal period"))
    {
        model->day_rash_onset[p->get_id()] = m->today();
        p->change_state(*m, RASH);
    }

    return ;

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_rash(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Checking if the agent will be detected or not
    bool detected = false;
    if (
        (m->par("Isolation period") >= 0) &&
        (m->runif() < 1.0/m->par("Days undetected"))
    )
    {
        model->agent_quarantine_triggered[p->get_id()] =
            QUARANTINE_PROCESS_ACTIVE;
        detected = true;
    }

    // Computing probabilities for state change
    m->array_double_tmp[0] = 1.0/m->par("Rash period"); // Recovery
    m->array_double_tmp[1] = m->par("Hospitalization rate"); // Hospitalization

    auto which = m->sample_from_probs(2);

    if (which == 0) // Recovers (probability 1/rash_period)
    {
        p->rm_virus(*m, detected ? ISOLATED_RECOVERED: RECOVERED);
    }
    else if (which == 1) // Hospitalized
    {
        m->record_hospitalization(*p);
        p->change_state(*m, detected ? DETECTED_HOSPITALIZED : HOSPITALIZED);
    }
    else if (which > 2)
    {
        throw std::logic_error("The roulette returned an unexpected value.");
    }
    else if (detected)
    {
        // If the agent is not hospitalized or recovered, then it is moved to
        // isolation.
        p->change_state(*m, ISOLATED);
        model->day_flagged[p->get_id()] = m->today();
    }

    return ;

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_isolated(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Figuring out if the agent can be released from isolation
    // if the isolation period is over.
    int days_since = m->today() - model->day_rash_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    // Sampling from the probabilities of recovery
    m->array_double_tmp[0] = 1.0/m->par("Rash period");

    // And hospitalization
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    auto which = m->sample_from_probs(2);

    // Recovers (which == 0 fires with probability 1/rash_period)
    if (which == 0)
    {
        p->rm_virus(*m, unisolate ? RECOVERED : ISOLATED_RECOVERED);
    }
    else if (which == 1)
    {

        m->record_hospitalization(*p);
        p->change_state(*m, unisolate ? HOSPITALIZED : DETECTED_HOSPITALIZED);

    }
    else if (unisolate)
    {
        p->change_state(*m, RASH);
    }


};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_quarantine_suscep(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Figuring out if the agent can be released from quarantine
    // if the quarantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    if (unquarantine)
    {
        p->change_state(*m, SUSCEPTIBLE);
    }

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_quarantine_latent(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Figuring out if the agent can be released from quarantine
    // if the quarantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    if (m->runif() < 1.0/(p->get_virus()->get_incubation(m)))
    {
        p->change_state(*m, unquarantine ? PRODROMAL : QUARANTINED_PRODROMAL);
    }
    else if (unquarantine)
    {
        p->change_state(*m, LATENT);
    }

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_quarantine_prodromal(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Otherwise, these are moved to the prodromal period, if
    // the quarantine period is over.
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

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_quarantine_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);
    int days_since = m->today() - model->day_flagged[p->get_id()];

    if (days_since >= m->par("Quarantine period"))
        p->change_state(*m, RECOVERED);

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_isolated_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto* model = model_cast<ModelMeaslesMixing<TSeq>,TSeq>(m);

    // Figuring out if the agent can be released from isolation
    // if the isolation period is over.
    int days_since = m->today() - model->day_rash_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    if (unisolate)
    {
        p->change_state(*m, RECOVERED);
    }

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_update_hospitalized(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // The agent is removed from the system
    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_virus(*m, RECOVERED);

};

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::m_quarantine_process() {

    // Process entity-level quarantine
    for (size_t agent_i = 0u; agent_i < this->size(); ++agent_i)
    {

        // Checking if the quarantine in the agent was triggered
        // or not
        if (agent_quarantine_triggered[agent_i] != QUARANTINE_PROCESS_ACTIVE)
            continue;

        if (Model<TSeq>::par("Quarantine period") < 0)
            continue;

        // Getting the number of contacts, if it is greater
        // than the maximum, it means that we overflowed, so
        // we will only quarantine the first EPI_MAX_TRACKING
        size_t n_contacts = contact_tracing.get_n_contacts(agent_i);
        if (n_contacts >= EPI_MAX_TRACKING)
            n_contacts = EPI_MAX_TRACKING;

        // When the rash onset started (this is for contact tracing)
        size_t day_rash_onset_agent_i = this->day_rash_onset[agent_i];

        for (size_t contact_i = 0u; contact_i < n_contacts; ++contact_i)
        {

            // Checking if the contact is within the contact tracing days prior
            auto [contact_id, contact_date] = contact_tracing.get_contact(agent_i, contact_i);
            bool within_days_prior =
                (day_rash_onset_agent_i - contact_date) <=
                this->par("Contact tracing days prior");
            if (!within_days_prior)
                continue;

            // Checking if we will detect the contact
            if (this->runif() > this->par("Contact tracing success rate"))
                continue;

            auto & agent = this->get_agent(contact_id);

            if (agent.get_state() > RASH)
                continue;

            // Agents with some tool won't be quarantined
            if (agent.get_n_tools() != 0u)
                continue;

            if (
                quarantine_willingness[contact_id] &&
                (this->par("Quarantine period") >= 0)
            )
            {

                switch (agent.get_state())
                {
                    case SUSCEPTIBLE:
                        agent.change_state(*this, QUARANTINED_SUSCEPTIBLE);
                        day_flagged[contact_id] = this->today();
                        break;
                    case LATENT:
                        agent.change_state(*this, QUARANTINED_LATENT);
                        day_flagged[contact_id] = this->today();
                        break;
                    case PRODROMAL:
                        agent.change_state(*this, QUARANTINED_PRODROMAL);
                        day_flagged[contact_id] = this->today();
                        break;
                    case RASH:
                        if (isolation_willingness[contact_id])
                        {
                            agent.change_state(*this, ISOLATED);
                            day_flagged[contact_id] = this->today();
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
        agent_quarantine_triggered[agent_i] = QUARANTINE_PROCESS_DONE;
    }

    return;
}

/**
 * @brief Template for a Measles model with population mixing, quarantine, and contact tracing
 *
 * @param model A ModelMeaslesMixing<TSeq> object where to set up the model.
 * @param n Number of agents in the population
 * @param prevalence Initial prevalence (proportion of infected individuals)
 * @param contact_rate Average number of contacts (interactions) per step
 * @param transmission_rate Probability of transmission per contact
 * @param vax_efficacy The efficacy of the vaccine
 * @param vax_reduction_recovery_rate The reduction in recovery rate due to the vaccine
 * @param incubation_period Average incubation period in days
 * @param prodromal_period Average prodromal period in days
 * @param rash_period Average rash period in days
 * @param contact_matrix Contact matrix specifying mixing patterns between population groups
 * @param hospitalization_rate Rate at which infected individuals are hospitalized
 * @param hospitalization_period Average duration of hospitalization in days
 * @param days_undetected Average number of days an infected individual remains undetected
 * @param quarantine_period Duration of quarantine in days for latent contacts
 * @param quarantine_willingness Proportion of individuals willing to comply with quarantine
 * @param isolation_willingness Proportion of individuals willing to self-isolate when detected
 * @param isolation_period Duration of isolation in days for detected infected individuals
 * @param prop_vaccinated Proportion of vaccinated agents
 * @param contact_tracing_success_rate Probability of successfully identifying contacts during tracing
 * @param contact_tracing_days_prior Number of days prior to detection for contact tracing
 */
template<typename TSeq>
inline ModelMeaslesMixing<TSeq>::ModelMeaslesMixing(
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double vax_efficacy,
    epiworld_double vax_reduction_recovery_rate,
    epiworld_double incubation_period,
    epiworld_double prodromal_period,
    epiworld_double rash_period,
    std::vector< double > contact_matrix,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double days_undetected,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_double isolation_willingness,
    epiworld_fast_int isolation_period,
    epiworld_double prop_vaccinated,
    epiworld_double contact_tracing_success_rate,
    epiworld_fast_uint contact_tracing_days_prior
    )
{

    // Setting up the contact matrix
    this->contact_matrix = contact_matrix;

    // Setting up parameters
    this->add_param(contact_rate, "Contact rate");
    this->add_param(transmission_rate, "Transmission rate");
    this->add_param(incubation_period, "Incubation period");
    this->add_param(prodromal_period, "Prodromal period");
    this->add_param(rash_period, "Rash period");
    this->add_param(hospitalization_rate, "Hospitalization rate");
    this->add_param(hospitalization_period, "Hospitalization period");
    this->add_param(days_undetected, "Days undetected");
    this->add_param(quarantine_period, "Quarantine period");
    this->add_param(quarantine_willingness, "Quarantine willingness");
    this->add_param(isolation_willingness, "Isolation willingness");
    this->add_param(isolation_period, "Isolation period");
    this->add_param(contact_tracing_success_rate, "Contact tracing success rate");
    this->add_param(contact_tracing_days_prior, "Contact tracing days prior");
    this->add_param(prop_vaccinated, "Vaccination rate");
    this->add_param(vax_efficacy, "Vax efficacy");
    this->add_param(vax_reduction_recovery_rate, "(IGNORED) Vax improved recovery");

    // state
    this->add_state("Susceptible", m_update_susceptible);
    this->add_state("Latent", m_update_latent);
    this->add_state("Prodromal", m_update_prodromal);
    this->add_state("Rash", m_update_rash);
    this->add_state("Isolated", m_update_isolated);
    this->add_state("Isolated Recovered", m_update_isolated_recovered);
    this->add_state("Detected Hospitalized", m_update_hospitalized);
    this->add_state("Quarantined Latent", m_update_quarantine_latent);
    this->add_state("Quarantined Susceptible", m_update_quarantine_suscep);
    this->add_state("Quarantined Prodromal", m_update_quarantine_prodromal);
    this->add_state("Quarantined Recovered", m_update_quarantine_recovered);
    this->add_state("Hospitalized", m_update_hospitalized);
    this->add_state("Recovered");

    // Global function
    this->queuing_off();

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus("Measles", prevalence, true);
    virus.set_state(LATENT, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Transmission rate");
    virus.set_prob_recovery("Rash period");
    virus.set_incubation("Incubation period");

    this->add_virus(virus);

    // Designing the vaccine
    ToolVaccine<TSeq> vaccine("MMR");

    vaccine.set_susceptibility_reduction(this->get_param("Vax efficacy"));

    vaccine.set_distribution(
        distribute_tool_randomly(prop_vaccinated, true)
    );

    this->add_tool(vaccine);

    this->queuing_off(); // No queuing need

    // Adding the empty population
    this->agents_empty_graph(n);

    this->set_name("Measles with Mixing and Quarantine");

}

template<typename TSeq>
inline ModelMeaslesMixing<TSeq> & ModelMeaslesMixing<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /* queue_ */
)
{

    Model<TSeq>::initial_states_fun =
        create_init_function_seir<TSeq>(proportions_)
        ;

    return *this;

}

template<typename TSeq>
inline void ModelMeaslesMixing<TSeq>::next()
{
    this->m_update_model();
    Model<TSeq>::next();
}

#undef MM
#endif
