#ifndef EPIWORLD_SIRD_H 
#define EPIWORLD_SIRD_H

#include "../model-bones.hpp"

/**
 * @brief Template for a Susceptible-Infected-Removed-Deceased (SIRD) model
 * 
 * ![Model Diagram](../assets/img/sird.png)
 * 
 * @ingroup death_compartmental
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRD : public Model<TSeq>
{
public:

    ModelSIRD() {};

    
    /**
     * @brief Constructs a new SIRD model with the given parameters.
     * 
     * @param vname The name of the vertex associated with this model.
     * @param prevalence The initial prevalence of the disease in the population.
     * @param transmission_rate The rate at which the disease spreads from infected to susceptible individuals.
     * @param recovery_rate The rate at which infected individuals recover and become immune.
     * @param death_rate The rate at which infected individuals die.
     */
    ///@{
    ModelSIRD(
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );
    ///@}

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with two elements:
     * - The proportion of non-infected individuals who have recovered.
     * - The proportion of non-infected individuals who have died.
    */
    ModelSIRD<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    );
    
};

template<typename TSeq>
inline ModelSIRD<TSeq>::ModelSIRD(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate, 
    epiworld_double death_rate
    )
{

    // Adding statuses
    this->add_state("Susceptible", default_update_susceptible<TSeq>);
    this->add_state("Infected", default_update_exposed<TSeq>);
    this->add_state("Recovered");
    this->add_state("Deceased");

    // Setting up parameters
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(transmission_rate, "Transmission rate"),
    this->add_param(death_rate, "Death rate");

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,2,3);
    virus.set_prob_recovery("Recovery rate");
    virus.set_prob_infecting("Transmission rate");
    virus.set_prob_death("Death rate");
    
    this->add_virus(virus);

    this->set_name("Susceptible-Infected-Recovered-Deceased (SIRD)");

   
}

template<typename TSeq>
inline ModelSIRD<TSeq> & ModelSIRD<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /**/ 
) {

    Model<TSeq>::initial_states_fun = 
        create_init_function_sird<TSeq>(proportions_)
        ;

    return *this;

}

#endif
