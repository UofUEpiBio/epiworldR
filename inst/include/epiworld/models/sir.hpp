#ifndef EPIWORLD_SIR_H 
#define EPIWORLD_SIR_H

#include "../model-bones.hpp"

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * ![Model Diagram](../assets/img/sir.png)
 * 
 * @ingroup basic_compartmental
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIR : public Model<TSeq>
{
public:

    ModelSIR() {};

    ModelSIR(
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate
    );

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with a single element:
     * - The proportion of non-infected individuals who have recovered.
    */
    ModelSIR<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    ) override;
    
};

template<typename TSeq>
inline ModelSIR<TSeq>::ModelSIR(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate
    )
{

    // Adding statuses
    this->add_state("Susceptible", default_update_susceptible<TSeq>);
    this->add_state("Infected", default_update_exposed<TSeq>);
    this->add_state("Recovered");

    // Setting up parameters
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(transmission_rate, "Transmission rate");

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,2,2);
    
    virus.set_prob_recovery("Recovery rate");
    virus.set_prob_infecting("Transmission rate");
    
    this->add_virus(virus);

    this->set_name("Susceptible-Infected-Recovered (SIR)");
   
}

template<typename TSeq>
inline ModelSIR<TSeq> & ModelSIR<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /* queue_ */
) {

    Model<TSeq>::initial_states_fun =
        create_init_function_sir<TSeq>(proportions_)
        ;

    return *this;

}

#endif
