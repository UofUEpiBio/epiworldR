#ifndef EPIWORLD_MODELS_SISD_HPP 
#define EPIWORLD_MODELS_SISD_HPP

/**
 * @brief Template for a Susceptible-Infected-Susceptible-Deceased (SISD) model
 * 
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 * @param inital_death epiworld_double Initial death_rate of the immune system
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSISD : public epiworld::Model<TSeq>
{

public:

    ModelSISD() {};

    ModelSISD(
        ModelSISD<TSeq> & model,
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );

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
    ModelSISD<TSeq> & model,
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    model.set_name("Susceptible-Infected-Susceptible-Deceased (SISD)");

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Infected", epiworld::default_update_exposed<TSeq>);
    model.add_state("Deceased");

    // Setting up parameters
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(recovery_rate, "Recovery rate");
    model.add_param(death_rate, "Death rate");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,0,2);
    
    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_prob_recovery(&model("Recovery rate"));
    virus.set_prob_death(0.01);
    
    model.add_virus(virus);

    return;

}

template<typename TSeq>
inline ModelSISD<TSeq>::ModelSISD(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    ModelSISD<TSeq>(
        *this,
        vname,
        prevalence,
        transmission_rate,
        recovery_rate,
        death_rate
    );    

    return;

}

#endif
