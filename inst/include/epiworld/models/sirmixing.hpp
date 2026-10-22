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
    public Model<TSeq>
{
private:

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

};

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

    // The contacts come from one mixing layer, grouped by entity
    auto mixing = layer::mixing<TSeq>(
        "mixing",
        std::vector< int >({ModelSIRMixing<TSeq>::INFECTED}),
        std::vector< int >({ModelSIRMixing<TSeq>::SUSCEPTIBLE, ModelSIRMixing<TSeq>::INFECTED, ModelSIRMixing<TSeq>::RECOVERED})
    );
    mixing->set_contact_matrix(contact_matrix);
    this->add_layer(std::move(mixing));

    UpdateFun<TSeq> update_susceptible = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void
        {

            Virus<TSeq> * virus = m->sample_virus(p);
            if (virus != nullptr)
                p->set_virus(*m, *virus, ModelSIRMixing<TSeq>::INFECTED);

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


    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(INFECTED, RECOVERED, RECOVERED);

    virus.set_prob_infecting("Prob. Transmission");
    virus.set_prob_recovery("Prob. Recovery");

    this->add_virus(virus);

    this->queuing_off(); // No queuing need

    // Adding the empty population
    this->set_population_size(n);

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
