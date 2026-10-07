#ifndef EPIWORLD_MODELS_SEIRMIXING_HPP
#define EPIWORLD_MODELS_SEIRMIXING_HPP

#include "../model-bones.hpp"

#define MM(i, j, n) \
    j * n + i

/**
 * @file seirentitiesconnected.hpp
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model with mixing
 * 
 * ![Model Diagram](../assets/img/seirmixing.png)
 * 
 * **Implementation details:**
 * <a href="../impl/mixing-and-entity-distribution.md">Mixing and Entity Distribution</a>,
 * <a href="../impl/sampling-contacts.md">Sampling Contacts</a>
 * 
 * @ingroup mixing_models
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRMixing :
    public Model<TSeq>,
    public ContactMatrix
{
private:

    // Samples the infectious contacts (group pools, available counts)
    sampler::Mixing<TSeq> mixing;

    #ifdef EPI_DEBUG
    std::vector< int > sampled_sizes;
    #endif

public:

    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int RECOVERED   = 3;

    ModelSEIRMixing() = delete;

    /**
     * @brief Constructs a ModelSEIRMixing object.
     *
     * @param vname The name of the ModelSEIRMixing object.
     * @param n The number of entities in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param avg_incubation_days The average incubation period of the disease in the model.
     * @param recovery_rate The recovery rate of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model. Each entry
     * (i,j) represents the expected number of contacts an agent in group i has
     * with agents in group j per day.
     */
    ModelSEIRMixing(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        std::vector< double > contact_matrix
    );

    void reset() override;

    std::unique_ptr< Model<TSeq> > clone_ptr() override;

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with a single element:
     * - The proportion of non-infected individuals who have recovered.
    */
    ModelSEIRMixing<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    ) override;

};

template<typename TSeq>
inline void ModelSEIRMixing<TSeq>::reset()
{

    Model<TSeq>::reset();

    // Checking contact matrix dimensions
    size_t nentities = this->entities.size();
    this->validate_contact_matrix(nentities);

    mixing.reset(*this);

    return;

}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSEIRMixing<TSeq>::clone_ptr()
{

    return std::make_unique<ModelSEIRMixing<TSeq>>(*this);

}


/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 *
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 * @param contact_matrix Contact matrix specifying expected contacts between groups.
 * Each entry (i,j) represents the expected number of contacts an agent in
 * group i has with agents in group j per day.
 */
template<typename TSeq>
inline ModelSEIRMixing<TSeq>::ModelSEIRMixing(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    std::vector< double > contact_matrix
    )
{

    mixing = sampler::Mixing<TSeq>(
        {ModelSEIRMixing<TSeq>::INFECTED},
        {ModelSEIRMixing<TSeq>::SUSCEPTIBLE, ModelSEIRMixing<TSeq>::EXPOSED, ModelSEIRMixing<TSeq>::INFECTED, ModelSEIRMixing<TSeq>::RECOVERED}
    );

    // Setting up the contact matrix
    this->set_contact_matrix(contact_matrix, true);

    UpdateFun<TSeq> update_susceptible = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void
        {

            if (p->get_n_entities() == 0)
                return;

            // Downcasting to retrieve the sampler attached to the
            // class
            auto * m_down = model_cast<ModelSEIRMixing<TSeq>, TSeq>(m);

            size_t ndraws = m_down->mixing.sample(p, *m, *m_down);

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

                auto & neighbor = m->get_agent(m_down->mixing.sampled_ids()[n]);

                auto & v = neighbor.get_virus();

                #ifdef EPI_DEBUG
                if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
                    throw std::logic_error(
                        "Trying to add an extra element to a temporal array outside of the range."
                    );
                #endif

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
                ModelSEIRMixing<TSeq>::EXPOSED
                );

            return;

        };

    UpdateFun<TSeq> update_exposed_and_infected = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSEIRMixing<TSeq>::EXPOSED)
            {

                // Getting the virus
                auto & v = p->get_virus();

                // Does the agent become infected?
                if (m->runif() < 1.0/(v->get_incubation(m)))
                {

                    p->change_state(*m, ModelSEIRMixing<TSeq>::INFECTED);
                    return;

                }


            } else if (state == ModelSEIRMixing<TSeq>::INFECTED)
            {


                // Odd: Die, Even: Recover
                epiworld_fast_uint n_events = 0u;
                auto & v = p->get_virus();

                // Recover
                m->array_double_tmp[n_events++] =
                    1.0 - (1.0 - v->get_prob_recovery(m)) *
                        (1.0 - p->get_recovery_enhancer(v, *m));

                #ifdef EPI_DEBUG
                if (n_events == 0u)
                {
                    printf_epiworld(
                        "[epi-debug] agent %i has 0 possible events!!\n",
                        static_cast<int>(p->get_id())
                        );
                    throw std::logic_error("Zero events in exposed.");
                }
                #else
                if (n_events == 0u)
                    return;
                #endif


                // Running the roulette
                int which = roulette(n_events, m);

                if (which < 0)
                    return;

                // Which roulette happen?
                p->rm_virus(*m);

                return ;

            } else
                throw std::logic_error("This function can only be applied to exposed or infected individuals. (SEIR)") ;

            return;

        };

    // Setting up parameters
    this->add_param(transmission_rate, "Prob. Transmission");
    this->add_param(recovery_rate, "Prob. Recovery");
    this->add_param(avg_incubation_days, "Avg. Incubation days");

    // state
    this->add_state("Susceptible", update_susceptible);
    this->add_state("Exposed", update_exposed_and_infected);
    this->add_state("Infected", update_exposed_and_infected);
    this->add_state("Recovered");

    // Global function
    GlobalFun<TSeq> update = [](Model<TSeq> * m) -> void
    {

        auto * m_down = model_cast<ModelSEIRMixing<TSeq>, TSeq>(m);

        m_down->mixing.update(*m);

        return;

    };

    this->add_globalevent(update, "Update infected individuals");


    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(EXPOSED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Prob. Transmission");
    virus.set_prob_recovery("Prob. Recovery");
    virus.set_incubation("Avg. Incubation days");

    this->add_virus(virus);

    this->queuing_off(); // No queuing need

    // Adding the empty population
    this->agents_empty_graph(n);

    this->set_name("Susceptible-Exposed-Infected-Removed (SEIR) with Mixing");

}

template<typename TSeq>
inline ModelSEIRMixing<TSeq> & ModelSEIRMixing<TSeq>::initial_states(
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
#endif
