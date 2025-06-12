#ifndef EPIWORLD_MODELS_SEIRMIXINGQUARANTINE_HPP
#define EPIWORLD_MODELS_SEIRMIXINGQUARANTINE_HPP

using namespace epiworld;

#define MM(i, j, n) \
    j * n + i

#if __cplusplus >= 202302L
    // C++23 or later
    #define GET_MODEL(model, output) \
        auto * output = dynamic_cast< ModelSEIRMixingQuarantine<TSeq> * >( (model) ); \
        /*Using the [[assume(...)]] to avoid the compiler warning \
        if the standard is C++23 or later */ \
        [[assume((output) != nullptr)]]
#else
    // C++17 or C++20
    #define GET_MODEL(model, output) \
        auto * output = dynamic_cast< ModelSEIRMixingQuarantine<TSeq> * >( (model) ); \
        assert((output) != nullptr); // Use assert for runtime checks
#endif

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
 * @file seirentitiesconnected.hpp
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model with mixing
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

    void m_update_infected_list();
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
    static void m_update_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_infected(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_suscep(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_quarantine_exposed(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_hospitalized(Agent<TSeq> * p, Model<TSeq> * m);
    static void m_update_isolated_recovered(Agent<TSeq> * p, Model<TSeq> * m);

    // Data about the quarantine process
    std::vector< bool > quarantine_willingness; ///< Indicator
    std::vector< bool > entity_quarantine_triggered; ///< Whether the quarantine process has started
    std::vector< bool > entity_can_quarantine; ///< Whether the entity can quarantine
    std::vector< int > day_flagged; ///< Either detected or started quarantine
    std::vector< int > day_onset; ///< Day of onset of the disease
    std::vector< int > day_exposed; ///< Day of exposure

    void m_quarantine_process();
    static void m_update_model(Model<TSeq> * m);

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

    ModelSEIRMixingQuarantine() {};
    
    /**
     * @brief Constructs a ModelSEIRMixingQuarantine object.
     *
     * @param model A reference to an existing ModelSEIRMixingQuarantine object.
     * @param vname The name of the ModelSEIRMixingQuarantine object.
     * @param n The number of entities in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param contact_rate The contact rate between entities in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param avg_incubation_days The average incubation period of the disease in the model.
     * @param recovery_rate The recovery rate of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model. Specified in
     * column-major order.
     */
    ModelSEIRMixingQuarantine(
        ModelSEIRMixingQuarantine<TSeq> & model,
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        std::vector< double > contact_matrix,
        std::vector< bool > entity_can_quarantine,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double days_undetected,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_fast_int isolation_period
    );
    
    /**
     * @brief Constructs a ModelSEIRMixingQuarantine object.
     *
     * @param vname The name of the ModelSEIRMixingQuarantine object.
     * @param n The number of entities in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param contact_rate The contact rate between entities in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param avg_incubation_days The average incubation period of the disease in the model.
     * @param recovery_rate The recovery rate of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model.
     */
    ModelSEIRMixingQuarantine(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        std::vector< double > contact_matrix,
        std::vector< bool > entity_can_quarantine,
        epiworld_double hospitalization_rate,
        epiworld_double hospitalization_period,
        // Policy parameters
        epiworld_double days_undetected,
        epiworld_fast_int quarantine_period,
        epiworld_double quarantine_willingness,
        epiworld_fast_int isolation_period
    );

    ModelSEIRMixingQuarantine<TSeq> & run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    void reset();

    Model<TSeq> * clone_ptr();

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with a single element:
     * - The proportion of non-infected individuals who have recovered.
    */
    ModelSEIRMixingQuarantine<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    );

    void set_contact_matrix(std::vector< double > cmat)
    {
        contact_matrix = cmat;
        return;
    };

    void set_entity_can_quarantine(std::vector< bool > can_quarantine)
    {
        entity_can_quarantine = can_quarantine;
        return;
    };

    std::vector< double > get_contact_matrix() const
    {
        return contact_matrix;
    };

    std::vector< bool > get_entity_can_quarantine() const
    {
        return entity_can_quarantine;
    };

    std::vector< bool > get_entity_quarantine_triggered() const
    {
        return entity_quarantine_triggered;
    };

    std::vector< bool > get_quarantine_willingness() const
    {
        return quarantine_willingness;
    };

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_model(Model<TSeq> * m)
{
    GET_MODEL(m, model);
    model->m_quarantine_process();
    model->events_run();
    model->m_update_infected_list();
    return;
}


template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_infected_list()
{

    auto & agents = Model<TSeq>::get_agents();

    std::fill(n_infected_per_group.begin(), n_infected_per_group.end(), 0u);
    
    for (auto & a : agents)
    {

        if (a.get_state() == ModelSEIRMixingQuarantine<TSeq>::INFECTED)
        {
            if (a.get_n_entities() > 0u)
            {
                const auto & entity = a.get_entity(0u);
                infected[
                    // Position of the group in the `infected` vector
                    entity_indices[entity.get_id()] +
                    // Position of the agent in the group
                    n_infected_per_group[entity.get_id()]++
                ] = a.get_id();

            }
        }

    }

    return;

}

template<typename TSeq>
inline size_t ModelSEIRMixingQuarantine<TSeq>::sample_agents(
    epiworld::Agent<TSeq> * agent,
    std::vector< size_t > & sampled_agents
    )
{

    size_t agent_group_id = agent->get_entity(0u).get_id();
    size_t ngroups = this->entities.size();

    int samp_id = 0;
    for (size_t g = 0; g < ngroups; ++g)
    {

        size_t group_size = n_infected_per_group[g];

        if (group_size == 0u)
            continue;

        // How many from this entity?
        int nsamples = epiworld::Model<TSeq>::rbinom(
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
            int which = epiworld::Model<TSeq>::runif() * group_size;

            // Correcting overflow error
            if (which >= static_cast<int>(group_size))
                which = static_cast<int>(group_size) - 1;

            #ifdef EPI_DEBUG
            auto & a = this->population.at(infected.at(entity_indices[g] + which));
            #else
            auto & a = this->get_agent(infected[entity_indices[g] + which]);
            #endif

            #ifdef EPI_DEBUG
            if (a.get_state() != ModelSEIRMixingQuarantine<TSeq>::INFECTED)
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
inline ModelSEIRMixingQuarantine<TSeq> & ModelSEIRMixingQuarantine<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    Model<TSeq>::run(ndays, seed);

    return *this;

}

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::reset()
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

    
    // Checking the quarantine variable
    if (entity_can_quarantine.size() != this->entities.size())
        throw std::length_error(
            std::string("The entity_can_quarantine vector must have the same size as the number of entities. ") +
            std::to_string(entity_can_quarantine.size()) +
            std::string(" != ") + std::to_string(this->entities.size()) +
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
    sampled_agents.resize(Model<TSeq>::size());

    // We only do it once
    n_infected_per_group.resize(this->entities.size(), 0u);
    std::fill(n_infected_per_group.begin(), n_infected_per_group.end(), 0u);

    // We are assuming one agent per entity
    infected.resize(Model<TSeq>::size());
    std::fill(infected.begin(), infected.end(), 0u);

    // This will say when do the groups start in the `infected` vector
    entity_indices.resize(this->entities.size(), 0u);
    std::fill(entity_indices.begin(), entity_indices.end(), 0u);
    for (size_t i = 1u; i < this->entities.size(); ++i)
    {

        entity_indices[i] +=
            this->entities[i - 1].size() +
            entity_indices[i - 1]
            ;

    }
    
    // Adjusting contact rate
    adjusted_contact_rate.clear();
    adjusted_contact_rate.resize(this->entities.size(), 0.0);

    for (size_t i = 0u; i < this->entities.size(); ++i)
    {
                
        adjusted_contact_rate[i] = 
            Model<TSeq>::get_param("Contact rate") /
                static_cast< epiworld_double > (this->get_entity(i).size());


        // Possibly correcting for a small number of agents
        if (adjusted_contact_rate[i] > 1.0)
            adjusted_contact_rate[i] = 1.0;

    }

    this->m_update_infected_list();

    // Setting up the quarantine parameters
    quarantine_willingness.resize(this->size(), false);
    for (size_t idx = 0; idx < quarantine_willingness.size(); ++idx)
        quarantine_willingness[idx] =
            Model<TSeq>::runif() < this->par("Quarantine willingness");

    entity_quarantine_triggered.resize(this->size(), false);
    std::fill(
        entity_quarantine_triggered.begin(),
        entity_quarantine_triggered.end(),
        false
    );

    day_flagged.resize(this->size(), 0);
    std::fill(
        day_flagged.begin(),
        day_flagged.end(),
        0
    );

    day_onset.resize(this->size(), 0);
    std::fill(
        day_onset.begin(),
        day_onset.end(),
        0
    );

    day_exposed.resize(this->size(), 0);
    std::fill(
        day_exposed.begin(),
        day_exposed.end(),
        0
    );
    
    return;

}

template<typename TSeq>
inline Model<TSeq> * ModelSEIRMixingQuarantine<TSeq>::clone_ptr()
{
    
    ModelSEIRMixingQuarantine<TSeq> * ptr = new ModelSEIRMixingQuarantine<TSeq>(
        *dynamic_cast<const ModelSEIRMixingQuarantine<TSeq>*>(this)
        );

    #if __cplusplus >= 202302L
        // C++23 or later
        [[assume(ptr != nullptr)]]
    #else
        // C++17 or C++20
        assert(ptr != nullptr); // Use assert for runtime checks
    #endif

    return dynamic_cast< Model<TSeq> *>(ptr);

}

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_susceptible(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    if (p->get_n_entities() == 0)
        return;

    // Downcasting to retrieve the sampler attached to the
    // class
    GET_MODEL(m, m_down);
    
    size_t ndraws = m_down->sample_agents(p, m_down->sampled_agents);

    #ifdef EPI_DEBUG
    m_down->sampled_sizes.push_back(static_cast<int>(ndraws));
    #endif

    if (ndraws == 0u)
        return;
    
    // Drawing from the set
    int nviruses_tmp = 0;
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
            
        /* And it is a function of susceptibility_reduction as well */ 
        m->array_double_tmp[nviruses_tmp] =
            (1.0 - p->get_susceptibility_reduction(v, m)) * 
            v->get_prob_infecting(m) * 
            (1.0 - neighbor.get_transmission_reduction(v, m)) 
            ; 
    
        m->array_virus_tmp[nviruses_tmp++] = &(*v);

    }

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return;

    p->set_virus(
        *m->array_virus_tmp[which],
        m,
        ModelSEIRMixingQuarantine<TSeq>::EXPOSED
        );

    return; 

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // Getting the virus
    auto & v = p->get_virus();

    // Does the agent become infected?
    if (m->runif() < 1.0/(v->get_incubation(m)))
    {

        p->change_state(m, ModelSEIRMixingQuarantine<TSeq>::INFECTED);

        GET_MODEL(m, model);
        model->day_onset[p->get_id()] = m->today();

        return;

    }

    return;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_infected(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    GET_MODEL(m, model);

    // Checking if the agent is detected
    bool detected = false;
    if (
        (m->par("Isolation period") >= 0) &&
        (m->runif() < 1.0/m->par("Days undetected")) &&
        (p->get_n_entities() != 0u) &&
        (model->entity_can_quarantine[p->get_entity(0u).get_id()]) 
    )
    {
        
        model->entity_quarantine_triggered[p->get_entity(0u).get_id()] = true;
        detected = true;

    }

    // Odd: Die, Even: Recover
    const auto & v = p->get_virus();
    m->array_double_tmp[0] = 1.0 - (1.0 - v->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(v, m));
    m->array_double_tmp[1] = m->par("Hospitalization rate");
        
    SAMPLE_FROM_PROBS(2, which);
    
    if (which == 0) // Recovers
    {
        if (detected)
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::ISOLATED_RECOVERED
            );
        }
        else
        {
            p->rm_virus(
                m, ModelSEIRMixingQuarantine<TSeq>::RECOVERED
            );
        }

        return;
    }
    else if (which == 1) // Hospitalized
    {

        if (detected)
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::DETECTED_HOSPITALIZED
            );
        }
        else
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::HOSPITALIZED
            );
        }

    }
    else if ((which == 2) && detected) // Nothing, but detected
    {
        // If the agent is detected, it goes to isolation
        p->change_state(
            m, ModelSEIRMixingQuarantine<TSeq>::ISOLATED
        );

    }
    
    return ;

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_isolated(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    GET_MODEL(m, model);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    // Sampling from the probabilities of recovery   
    m->array_double_tmp[0] = 1.0 -
        (1.0 - p->get_virus()->get_prob_recovery(m)) *
        (1.0 - p->get_recovery_enhancer(p->get_virus(), m));

    // And hospitalization
    m->array_double_tmp[1] = m->par("Hospitalization rate");

    SAMPLE_FROM_PROBS(2, which);

    // Recovers
    if (which == 0)
    {
        if (unisolate)
        {
            p->rm_virus(
                m,
                ModelSEIRMixingQuarantine<TSeq>::RECOVERED
            );
        }
        else
            p->rm_virus(
                m, ModelSEIRMixingQuarantine<TSeq>::ISOLATED_RECOVERED
            );
    }
    else if (which == 1)
    {

        if (unisolate)
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::HOSPITALIZED
            );
        }
        else
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::DETECTED_HOSPITALIZED
            );
        }
    }
    else if ((which == 2) && unisolate)
    {
        p->change_state(
            m, ModelSEIRMixingQuarantine<TSeq>::INFECTED
        );
    }


};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_quarantine_suscep(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    GET_MODEL(m, model);

    // Figuring out if the agent can be released from quarantine
    // if the quarantine period is over.
    int days_since = m->today() - model->day_flagged[p->get_id()];

    bool unquarantine =
        (m->par("Quarantine period") <= days_since) ?
        true: false;

    if (unquarantine)
    {
        p->change_state(
            m, ModelSEIRMixingQuarantine<TSeq>::SUSCEPTIBLE
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_quarantine_exposed(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    GET_MODEL(m, model);

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
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::INFECTED
            );
        }
        else
        {
            p->change_state(
                m, ModelSEIRMixingQuarantine<TSeq>::ISOLATED
            );
        }

    }
    else if (unquarantine)
    {
        p->change_state(
            m, ModelSEIRMixingQuarantine<TSeq>::EXPOSED
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_isolated_recovered(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    GET_MODEL(m, model);

    // Figuring out if the agent can be released from isolation
    // if the quarantine period is over.
    int days_since = m->today() - model->day_onset[p->get_id()];

    bool unisolate =
        (m->par("Isolation period") <= days_since) ?
        true: false;

    if (unisolate)
    {
        p->change_state(
            m, ModelSEIRMixingQuarantine<TSeq>::RECOVERED
        );
    }

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_update_hospitalized(
    Agent<TSeq> * p, Model<TSeq> * m
) {

    // The agent is removed from the system
    if (m->runif() < 1.0/m->par("Hospitalization period"))
        p->rm_virus(m, ModelSEIRMixingQuarantine<TSeq>::RECOVERED);

};

template<typename TSeq>
inline void ModelSEIRMixingQuarantine<TSeq>::m_quarantine_process() {

    int entity_num = -1;
    for (auto & entity: Model<TSeq>::get_entities())
    {

        // Checking if the quarantine in the entity was triggered
        // or not
        if (!entity_quarantine_triggered[++entity_num])
            return;

        if (
            (this->par("Quarantine period") < 0) &&
            (this->par("Isolation period") < 0)
        )
            return;

        for (auto i: entity.get_agents())
        {

            auto & agent = Model<TSeq>::get_agent(i);

            if (agent.get_state() > INFECTED)
                continue;

            // Agents with some tool won't be quarantined
            if (agent.get_n_tools() != 0u)
                continue;

            if (
                quarantine_willingness[i] &&
                (Model<TSeq>::par("Quarantine period") >= 0))
            {

                switch (agent.get_state())
                {
                    case SUSCEPTIBLE:
                        agent.change_state(this, QUARANTINED_SUSCEPTIBLE);
                        break;
                    case EXPOSED:
                        agent.change_state(this, QUARANTINED_EXPOSED);
                        break;
                    case INFECTED:
                        agent.change_state(this, ISOLATED);
                        break;
                    default:
                        throw std::logic_error(
                            "The agent is not in a state that can be quarantined."
                        );
                }

                // And we add the day of quarantine
                day_flagged[i] = Model<TSeq>::today();

            }
        }

        // Setting the quarantine process off
        entity_quarantine_triggered[entity_num] = false;
    }

    return;
}

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 */
template<typename TSeq>
inline ModelSEIRMixingQuarantine<TSeq>::ModelSEIRMixingQuarantine(
    ModelSEIRMixingQuarantine<TSeq> & model,
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    std::vector< double > contact_matrix,
    std::vector< bool > entity_can_quarantine,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double days_undetected,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_fast_int isolation_period
    )
{

    // Setting up the contact matrix
    this->contact_matrix = contact_matrix;
    this->entity_can_quarantine = entity_can_quarantine;

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Prob. Transmission");
    model.add_param(recovery_rate, "Prob. Recovery");
    model.add_param(avg_incubation_days, "Avg. Incubation days");
    model.add_param(hospitalization_rate, "Hospitalization rate");
    model.add_param(hospitalization_period, "Hospitalization period");
    model.add_param(days_undetected, "Days undetected");
    model.add_param(quarantine_period, "Quarantine period");
    model.add_param(
        quarantine_willingness, "Quarantine willingness"
    );
    model.add_param(isolation_period, "Isolation period");
    
    // state
    model.add_state("Susceptible", m_update_susceptible);
    model.add_state("Exposed", m_update_exposed);
    model.add_state("Infected", m_update_infected);
    model.add_state("Isolated", m_update_isolated);
    model.add_state("Detected Hospitalized", m_update_hospitalized);
    model.add_state("Quarantined Susceptible", m_update_quarantine_suscep);
    model.add_state("Quarantined Exposed", m_update_quarantine_exposed);
    model.add_state("Isolated Recovered", m_update_isolated_recovered);
    model.add_state("Hospitalized", m_update_hospitalized);
    model.add_state("Recovered");

    // Global function
    model.add_globalevent(this->m_update_model, "Update infected individuals");
    model.queuing_off();

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(
        ModelSEIRMixingQuarantine<TSeq>::EXPOSED,
        ModelSEIRMixingQuarantine<TSeq>::RECOVERED,
        ModelSEIRMixingQuarantine<TSeq>::RECOVERED
        );

    virus.set_prob_infecting(&model("Prob. Transmission"));
    virus.set_prob_recovery(&model("Prob. Recovery"));
    virus.set_incubation(&model("Avg. Incubation days"));

    model.add_virus(virus);

    model.queuing_off(); // No queuing need

    // Adding the empty population
    model.agents_empty_graph(n);

    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR) with Mixing");

    return;

}

template<typename TSeq>
inline ModelSEIRMixingQuarantine<TSeq>::ModelSEIRMixingQuarantine(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    std::vector< double > contact_matrix,
    std::vector< bool > entity_can_quarantine,
    epiworld_double hospitalization_rate,
    epiworld_double hospitalization_period,
    // Policy parameters
    epiworld_double days_undetected,
    epiworld_fast_int quarantine_period,
    epiworld_double quarantine_willingness,
    epiworld_fast_int isolation_period
    )
{   

    this->contact_matrix = contact_matrix;
    this->entity_can_quarantine = entity_can_quarantine;

    ModelSEIRMixingQuarantine(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        transmission_rate,
        avg_incubation_days,
        recovery_rate,
        contact_matrix,
        entity_can_quarantine,
        hospitalization_rate,
        hospitalization_period,
        // Policy parameters
        days_undetected,
        quarantine_period,
        quarantine_willingness,
        isolation_period
    );

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
#undef MM
#undef GET_MODEL
#endif
