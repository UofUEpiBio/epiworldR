#ifndef EPIWORLD_MODELS_SISD_HPP 
#define EPIWORLD_MODELS_SISD_HPP

#include "../model-bones.hpp"

/**
 * @brief Template for a Susceptible-Infected-Susceptible-Deceased (SISD) model
 * 
 * ![Model Diagram](../assets/img/sisd.png)
 * 
 * @ingroup death_compartmental
 * 
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 * @param inital_death epiworld_double Initial death_rate of the immune system
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSISD : public Model<TSeq>
{

public:

    ModelSISD() = delete;

    ModelSISD(
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );

};

template<typename TSeq>
inline ModelSISD<TSeq>::ModelSISD(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    this->set_name("Susceptible-Infected-Susceptible-Deceased (SISD)");

    // Adding statuses
    this->add_state("Susceptible", default_update_susceptible<TSeq>);
    this->add_state("Infected", default_update_exposed<TSeq>);
    this->add_state("Deceased");

    // Setting up parameters
    this->add_param(transmission_rate, "Transmission rate");
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(death_rate, "Death rate");

    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,0,2);
    
    virus.set_prob_infecting("Transmission rate");
    virus.set_prob_recovery("Recovery rate");
    virus.set_prob_death(0.01);
    
    this->add_virus(virus);

}

#endif
