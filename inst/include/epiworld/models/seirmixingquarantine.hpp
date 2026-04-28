#ifndef EPIWORLD_MODELS_SEIRMIXINGQUARANTINE_HPP
#define EPIWORLD_MODELS_SEIRMIXINGQUARANTINE_HPP

#include "../model-bones.hpp"

#define MM(i, j, n) \
    j * n + i

/**
 * @file seirmixingquarantine.hpp
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model with mixing, quarantine, and contact tracing
 */

/**
 * @brief SEIR model with mixing, quarantine, and contact tracing
 *
 * This class implements a Susceptible-Exposed-Infected-Removed (SEIR) epidemiological model
 * with additional features including:
 *
 * - Population mixing based on contact matrices
 * - Quarantine measures for exposed contacts
 * - Isolation policies for detected infected individuals
 * - Contact tracing with configurable success rates
 * - Hospitalization of severe cases
 * - Individual willingness to comply with public health measures
 *
 * The model supports 9 distinct states:
 *
 * - Susceptible: Individuals who can become infected
 * - Exposed: Infected but not yet infectious (incubation period)
 * - Infected: Infectious individuals in the community
 * - Isolated: Detected infected individuals in self-isolation
 * - Quarantined Susceptible: Susceptible individuals in quarantine due to contact tracing
 * - Quarantined Exposed: Exposed individuals in quarantine due to contact tracing
 * - Isolated Recovered: Recovered individuals still in isolation
 * - Hospitalized: Individuals requiring hospital care
 * - Recovered: Individuals who have recovered and gained immunity
 *
 * ![Model Diagram](../assets/img/seirmixingquarantine.png)
 *
 * @tparam TSeq Type for genetic sequences (default: EPI_DEFAULT_TSEQ)
 * @ingroup mixing_models
 * 
 * **Implementation details:**
 * <a href="../impl/quarantine-isolation-and-contact-tracing.md">Quarantine, Isolation, and Contact Tracing</a>,
 * <a href="../impl/mixing-and-entity-distribution.md">Mixing and Entity Distribution</a>,
 * <a href="../impl/sampling-contacts.md">Sampling Contacts</a>
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRMixingQuarantine : public Model<TSeq>
{
private:

    // Vector of vectors of infected agents
    std::vector< size_t > infected;

    // Number of infected agents in each group
    std::vector< size_t > n_infected_per_group;

    // Where the agents start in the `infected` vector
    std::vector< size_t > entity_indices;

    void _update_infected_list();
    std::vector< size_t > sampled_agents;
    size_t _sample_agents(
        Agent<TSeq> * agent,
        std::vector< size_t > & sampled_agents
        );
    std::vector< double > adjusted_contact_rate;
    std::vector< double > contact_matrix;

    #ifdef EPI_DEBUG
    std::vector< int > sampled_sizes;
    #endif

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
    std::vector< int > day_exposed; ///< Day of exposure

    static void _quarantine_process(Model<TSeq> * m);

public:

    static const int SUSCEPTIBLE             = 0;
    static const int EXPOSED                 = 1;
    static const int INFECTED                = 2;
    static const int ISOLATED                = 3;
    static const int QUARANTINED_SUSCEPTIBLE = 4;
    static const int QUARANTINED_EXPOSED     = 5;
    static const int ISOLATED_RECOVERED      = 6;
    static const int HOSPITALIZED            = 7;
    static const int RECOVERED               = 8;

    static const size_t QUARANTINE_PROCESS_INACTIVE = 0u;
    static const size_t QUARANTINE_PROCESS_ACTIVE   = 1u;
    static const size_t QUARANTINE_PROCESS_DONE     = 2u;

    ModelSEIRMixingQuarantine() = delete;

    /**
     * @brief Constructs a ModelSEIRMixingQuarantine object.
     *
     * @param vname The name of the ModelSEIRMixingQuarantine object.
     * @param n The number of agents in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param avg_incubation_days The average incubation period of the disease in the model.
     * @param recovery_rate The recovery rate of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model.
     * Specified in column-major order. Each entry (i,j) represents the
     * expected number of contacts an agent in group i has with agents in
     * group j per day.
     * @param hospitalization_rate The rate at which infected individuals are hospitalized.
     * @param hospitalization_period The average duration of hospitalization in days.
     * @param days_undetected The average number of days an infected individual remains undetected.
     * @param quarantine_period The duration of quarantine in days for exposed contacts.
     * @param quarantine_willingness The proportion of individuals willing to comply with quarantine measures.
     * @param isolation_willingness The proportion of individuals willing to self-isolate when detected.
     * @param isolation_period The duration of isolation in days for detected infected individuals.
     * @param contact_tracing_success_rate The probability of successfully identifying and tracing contacts (default: 1.0).
     * @param contact_tracing_days_prior The number of days prior to detection for which contacts are traced (default: 4).
     */
    ModelSEIRMixingQuarantine(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        std::vector< double > contact_matrix,
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
    ModelSEIRMixingQuarantine<TSeq> & initial_states(
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

    void next() override;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_infected_list()
{

    auto & agents = this->get_agents();

    std::fill(n_infected_per_group.begin(), n_infected_per_group.end(), 0u);

    // Resetting the number of available contacts
    adjusted_contact_rate.assign(this->entities.size(), 0.0);

    for (const auto & a : agents)
    {

        if (a.get_state() == INFECTED)
        {
            if (a.get_n_entities() > 0u)
            {
                const auto & entity = a.get_entity(0u, *this);
                infected[
                    // Position of the group in the `infected` vector
                    entity_indices[entity.get_id()] +
                    // Position of the agent in the group
                    n_infected_per_group[entity.get_id()]++
                ] = a.get_id();

            }
        }

        // Setting how many agents are available for contact
        if (
            ((a.get_state() < ISOLATED) || (a.get_state() == RECOVERED)) &&
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
            rate = 1.0 / rate;
        else
            rate = 0.0;  // No available contacts in this group

        if (rate > 1.0)
            rate = 1.0;
    }

    return;

}

template<typename TSeq>
inline size_t ModelSEIRMixingQuarantine<TSeq>::_sample_agents(
    Agent<TSeq> * agent,
    std::vector< size_t > & sampled_agents
    )
{

    size_t agent_group_id = agent->get_entity(0u, *this).get_id();
    size_t ngroups = this->entities.size();

    int samp_id = 0;
    for (size_t g = 0; g < ngroups; ++g)
    {

        size_t group_size = n_infected_per_group[g];

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
            auto & a = this->population.at(infected.at(entity_indices[g] + which));
            #else
            auto & a = this->get_agent(infected[entity_indices[g] + which]);
            #endif

            #ifdef EPI_DEBUG
            if (a.get_state() != INFECTED)
                throw std::logic_error(
                    "The agent is not infected, but it should be."
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
inline void ModelSEIRMixingQuarantine<TSeq>::reset()
{

    Model<TSeq>::reset();

    // Checking contact matrix dimensions
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
        for (size_t j = 0u; j < this->entities.size(); ++j)
        {
            if (this->contact_matrix[MM(i, j, nentities)] < 0.0)
                throw std::range_error(
                    std::string("The contact matrix must be non-negative. ") +
                    std::to_string(this->contact_matrix[MM(i, j, nentities)]) +
                    std::string(" < 0.")
                    );
        }
    }

    // Do it the first time only
    sampled_agents.resize(this->size());

    // We only do it once
    n_infected_per_group.assign(this->entities.size(), 0u);

    // We are assuming one agent per entity
    infected.assign(this->size(), 0u);

    // This will say when do the groups start in the `infected` vector
    entity_indices.assign(this->entities.size(), 0u);
    for (size_t i = 1u; i < this->entities.size(); ++i)
    {

        entity_indices[i] +=
            this->entities[i - 1].size() +
            entity_indices[i - 1]
            ;

    }

    this->_update_infected_list();

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
    day_exposed.assign(this->size(), 0);

    return;

}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSEIRMixingQuarantine<TSeq>::clone_ptr()
{

    return std::make_unique<ModelSEIRMixingQuarantine<TSeq>>(*this);

}

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_susceptible(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    if (p->get_n_entities() == 0)
        return;

    // Downcasting to retrieve the sampler attached to the
    // class
    auto * m_down = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    size_t ndraws = m_down->_sample_agents(p, m_down->sampled_agents);

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
        m_down->get_contact_tracing().add_contact(neighbor.get_id(), p->get_id(), m->today());

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

    p->set_virus(*m, 
        *m->array_virus_tmp[which],
        ModelSEIRMixingQuarantine<TSeq>::EXPOSED
        );

    return;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // Getting the virus
    auto & v = p->get_virus();

    // Does the agent become infected?
    if (m->runif() < 1.0/(v->get_incubation(m)))
    {

        p->change_state(*m, INFECTED);

        auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);
        model->day_onset[p->get_id()] = m->today();

    }

    return;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_infected(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Sampling whether the agent is detected or not.
    // If Days undetected < 0, detection is disabled (never detected).
    // If Days undetected == 0, the agent is always detected.
    epiworld_double days_undetected = m->par("Days undetected");
    bool detected = (days_undetected < 0.0) ?
        false : ((days_undetected == 0.0) ?
            true : (m->runif() < 1.0 / days_undetected));

    // If detected and the entity can quarantine, we start
    // the quarantine process
    if (detected)
    {
        model->agent_quarantine_triggered[p->get_id()] =
            QUARANTINE_PROCESS_ACTIVE;
    }

    // Checking if the agent is willing to isolate individually
    // This is separate from quarantine and can happen even if agent cannot quarantine
    bool isolation_detected = (m->par("Isolation period") >= 0) &&
        detected &&
        (model->isolation_willingness[p->get_id()])
    ;

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
            p->change_state(*m, ISOLATED_RECOVERED);
        }
        else
        {
            p->rm_virus(*m, RECOVERED);
        }

        return;
    }
    else if (which == 1) // Hospitalized
    {
        m->record_hospitalization(*p);
        p->change_state(*m, HOSPITALIZED);

    }
    else if ((which == 2) && isolation_detected) // Nothing, but detected
    {
        // If the agent is detected, it goes to isolation
        p->change_state(*m, ISOLATED);

    }

    return ;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_isolated(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

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
            p->rm_virus(*m, RECOVERED);
        }
        else
            p->rm_virus(*m, ISOLATED_RECOVERED);
    }
    else if (which == 1)
    {
        m->record_hospitalization(*p);
        p->change_state(*m, HOSPITALIZED);
    }
    else if ((which == 2) && unisolate)
    {
        p->change_state(*m, INFECTED);
    }


};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_quarantine_suscep(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Figuring out if the agent can be released from quarantine
    // if the quarantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    if (unquarantine)
    {
        p->change_state(*m, 
            ModelSEIRMixingQuarantine<TSeq>::SUSCEPTIBLE
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_quarantine_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Figuring out if the agent can be released from quarantine
    // if the quarantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    if (m->runif() < 1.0/(p->get_virus()->get_incubation(m)))
    {

        // Recording the day of onset
        model->day_onset[p->get_id()] = m->today();

        // If the agent is unquarantined, it becomes infected
        if (unquarantine)
        {
            p->change_state(*m, 
                ModelSEIRMixingQuarantine<TSeq>::INFECTED
            );
        }
        else
        {
            p->change_state(*m, 
                ModelSEIRMixingQuarantine<TSeq>::ISOLATED
            );
        }

    }
    else if (unquarantine)
    {
        p->change_state(*m, 
            ModelSEIRMixingQuarantine<TSeq>::EXPOSED
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_isolated_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    if (unisolate)
    {
        p->change_state(*m, 
            ModelSEIRMixingQuarantine<TSeq>::RECOVERED
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_update_hospitalized(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // The agent is removed from the system
    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_virus(*m, ModelSEIRMixingQuarantine<TSeq>::RECOVERED);

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::_quarantine_process(Model<TSeq> * m) {

    auto * model = model_cast<ModelSEIRMixingQuarantine<TSeq>, TSeq>(m);

    // Process entity-level quarantine
    for (size_t agent_i = 0u; agent_i < m->size(); ++agent_i)
    {

        // Checking if the quarantine in the agent was triggered
        // or not
        if (
            model->agent_quarantine_triggered[agent_i] !=
            ModelSEIRMixingQuarantine<TSeq>::QUARANTINE_PROCESS_ACTIVE
        )
            continue;

        if (m->par("Quarantine period") < 0)
        {
            model->agent_quarantine_triggered[agent_i] =
            ModelSEIRMixingQuarantine<TSeq>::QUARANTINE_PROCESS_DONE;
            continue;
        }

        // Getting the number of contacts, if it is greater
        // than the maximum, it means that we overflowed, so
        // we will only quarantine the first EPI_MAX_TRACKING
        auto & ct = m->get_contact_tracing();
        size_t n_contacts = ct.get_n_contacts(agent_i);
        if (n_contacts >= EPI_MAX_TRACKING)
            n_contacts = EPI_MAX_TRACKING;

        auto success_rate = m->par("Contact tracing success rate");
        auto days_prior = m->par("Contact tracing days prior");
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
                        agent.change_state(*m, ModelSEIRMixingQuarantine<TSeq>::QUARANTINED_SUSCEPTIBLE);
                        model->day_flagged[contact_id] = m->today();
                        break;
                    case EXPOSED:
                        agent.change_state(*m, ModelSEIRMixingQuarantine<TSeq>::QUARANTINED_EXPOSED);
                        model->day_flagged[contact_id] = m->today();
                        break;
                    case INFECTED:
                        if (model->isolation_willingness[contact_id])
                        {
                            agent.change_state(*m, ModelSEIRMixingQuarantine<TSeq>::ISOLATED);
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
        model->agent_quarantine_triggered[agent_i] =
            ModelSEIRMixingQuarantine<TSeq>::QUARANTINE_PROCESS_DONE;
    }

    return;
}

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model with mixing, quarantine, and contact tracing
 *
 * @param model A ModelSEIRMixingQuarantine<TSeq> object where to set up the SEIR model.
 * @param vname Name of the virus
 * @param n Number of agents in the population
 * @param prevalence Initial prevalence (proportion of infected individuals)
 * @param transmission_rate Probability of transmission per contact
 * @param avg_incubation_days Average incubation period in days
 * @param recovery_rate Probability of recovery per day
 * @param contact_matrix Contact matrix specifying expected contacts between groups.
 * Each entry (i,j) represents the expected number of contacts an agent in
 * group i has with agents in group j per day.
 * @param hospitalization_rate Rate at which infected individuals are hospitalized
 * @param hospitalization_period Average duration of hospitalization in days
 * @param days_undetected Average number of days an infected individual remains undetected
 * @param quarantine_period Duration of quarantine in days for exposed contacts
 * @param quarantine_willingness Proportion of individuals willing to comply with quarantine
 * @param isolation_willingness Proportion of individuals willing to self-isolate when detected
 * @param isolation_period Duration of isolation in days for detected infected individuals
 * @param contact_tracing_success_rate Probability of successfully identifying contacts during tracing
 * @param contact_tracing_days_prior Number of days prior to detection for contact tracing
 */
template<typename TSeq>
inline ModelSEIRMixingQuarantine<TSeq>::ModelSEIRMixingQuarantine(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    std::vector< double > contact_matrix,
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

    // Setting up the contact matrix
    this->contact_matrix = contact_matrix;

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

    // state
    this->add_state("Susceptible", _update_susceptible);
    this->add_state("Exposed", _update_exposed);
    this->add_state("Infected", _update_infected);
    this->add_state("Isolated", _update_isolated);
    this->add_state("Quarantined Susceptible", _update_quarantine_suscep);
    this->add_state("Quarantined Exposed", _update_quarantine_exposed);
    this->add_state("Isolated Recovered", _update_isolated_recovered);
    this->add_state("Hospitalized", _update_hospitalized);
    this->add_state("Recovered");

    // Global function
    this->add_globalevent(
        this->_quarantine_process,
        "Update infected individuals"
    );

    this->queuing_off();

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(EXPOSED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Prob. Transmission");
    virus.set_prob_recovery("Prob. Recovery");
    virus.set_incubation("Avg. Incubation days");

    this->add_virus(virus);

    this->queuing_off(); // No queuing need

    // Enable contact tracing for quarantine process
    this->contact_tracing_on(EPI_MAX_TRACKING);

    // Adding the empty population
    this->agents_empty_graph(n);

    this->set_name("SEIR with Mixing and Quarantine");

    return;

}

template<typename TSeq>
inline ModelSEIRMixingQuarantine<TSeq> & ModelSEIRMixingQuarantine<TSeq>::initial_states(
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
inline void ModelSEIRMixingQuarantine<TSeq>::next() {

    this->_update_infected_list();
    Model<TSeq>::next();

}

#undef MM
#endif
