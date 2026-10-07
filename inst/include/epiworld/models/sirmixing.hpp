#ifndef EPIWORLD_MODELS_SIRMIXING_HPP
#define EPIWORLD_MODELS_SIRMIXING_HPP

#include "../model-bones.hpp"

/**
 * @file sirmixing.hpp
 * @brief Template for a Susceptible-Infected-Removed (SIR) model with mixing
 * 
 * ![Model Diagram](../assets/img/sirmixing.png)
 * 
 * **Implementation details:**
 * <a href="../impl/mixing-and-entity-distribution.md">Mixing and Entity Distribution</a>,
 * <a href="../impl/sampling-contacts.md">Sampling Contacts</a>
 * 
 * @ingroup mixing_models
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRMixing :
    public Model<TSeq>,
    public ContactMatrix
{
private:

    // Samples the infectious contacts (group pools, available counts)
    sampler::Mixing<TSeq> mixing;

    size_t index(size_t i, size_t j, size_t n) {
        return j * n + i;
    }

public:

    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;

    ModelSIRMixing() = delete;

    /**
     * @brief Constructs a ModelSIRMixing object.
     *
     * @param vname The name of the ModelSIRMixing object.
     * @param n The number of agents in the model.
     * @param prevalence The initial prevalence of the disease in the model.
     * @param transmission_rate The transmission rate of the disease in the model.
     * @param recovery_rate The recovery rate of the disease in the model.
     * @param contact_matrix The contact matrix between entities in the model.
     * Specified in column-major order. Each entry (i,j) represents the
     * expected number of contacts an agent in group i has with agents in
     * group j per day. Entry (i,j) is located at `contact_matrix[j * n + i]`.
     */
    ModelSIRMixing(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
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
    ModelSIRMixing<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    ) override;

    size_t get_n_infected(size_t group) const
    {
        return mixing.get_n_infectious(group);
    }

};

template<typename TSeq>
inline void ModelSIRMixing<TSeq>::reset()
{

    Model<TSeq>::reset();

    // Checking contact matrix dimensions
    size_t nentities = this->entities.size();
    this->validate_contact_matrix(nentities);

    mixing.reset(*this);

    return;
}

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSIRMixing<TSeq>::clone_ptr()
{

    return std::make_unique<ModelSIRMixing<TSeq>>(*this);

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
inline ModelSIRMixing<TSeq>::ModelSIRMixing(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    std::vector< double > contact_matrix
    )
{

    mixing = sampler::Mixing<TSeq>(
        {ModelSIRMixing<TSeq>::INFECTED},
        {ModelSIRMixing<TSeq>::SUSCEPTIBLE, ModelSIRMixing<TSeq>::INFECTED, ModelSIRMixing<TSeq>::RECOVERED}
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
            auto * m_down = model_cast<ModelSIRMixing<TSeq>, TSeq>(m);

            size_t ndraws = m_down->mixing.sample(p, *m, *m_down);

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
                    throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
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
                ModelSIRMixing<TSeq>::INFECTED
                );

            return;

        };

    UpdateFun<TSeq> update_infected = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSIRMixing<TSeq>::INFECTED)
            {


                // Odd: Die, Even: Recover
                epiworld_fast_uint n_events = 0u;
                auto & v = p->get_virus();

                // Recover
                m->array_double_tmp[n_events++] =
                    1.0 - (1.0 - v->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(v, *m));

                #ifdef EPI_DEBUG
                if (n_events == 0u)
                {
                    printf_epiworld(
                        "[epi-debug] agent %i has 0 possible events!!\n",
                        static_cast<int>(p->get_id())
                        );
                    throw std::logic_error("Zero events in infected.");
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
                throw std::logic_error("This function can only be applied to infected individuals. (SIR)") ;

            return;

        };

    // Setting up parameters
    this->add_param(transmission_rate, "Prob. Transmission");
    this->add_param(recovery_rate, "Prob. Recovery");

    // state
    this->add_state("Susceptible", update_susceptible);
    this->add_state("Infected", update_infected);
    this->add_state("Recovered");

    // Global function
    GlobalFun<TSeq> update = [](Model<TSeq> * m) -> void
    {

        auto * m_down = model_cast<ModelSIRMixing<TSeq>, TSeq>(m);

        m_down->mixing.update(*m);

        return;

    };

    this->add_globalevent(update, "Update infected individuals");


    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(INFECTED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Prob. Transmission");
    virus.set_prob_recovery("Prob. Recovery");

    this->add_virus(virus);

    this->queuing_off(); // No queuing need

    // Adding the empty population
    this->agents_empty_graph(n);

    this->set_name("Susceptible-Infected-Removed (SIR) with Mixing");

    return;

}

template<typename TSeq>
inline ModelSIRMixing<TSeq> & ModelSIRMixing<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /* queue_ */
)
{

    Model<TSeq>::initial_states_fun =
        create_init_function_sir<TSeq>(proportions_)
        ;

    return *this;

}

#endif
