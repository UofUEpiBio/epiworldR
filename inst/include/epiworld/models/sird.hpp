#ifndef EPIWORLD_SIRD_H 
#define EPIWORLD_SIRD_H

/**
 * @brief Template for a Susceptible-Infected-Removed-Deceased (SIRD) model
 * 
 * @param model A Model<TSeq> object where to set up the SIRD.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 * @param initial_death epiworld_double Initial death_rate of the immune system 
 */
template<typename TSeq = int>
class ModelSIRD : public epiworld::Model<TSeq>
{
public:

    ModelSIRD() {};

    ModelSIRD(
        ModelSIRD<TSeq> & model,
        std::string vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );

    ModelSIRD(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );
    
};

template<typename TSeq>
inline ModelSIRD<TSeq>::ModelSIRD(
    ModelSIRD<TSeq> & model,
    std::string vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate, 
    epiworld_double death_rate
    )
{

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Infected", epiworld::default_update_exposed<TSeq>);
    model.add_state("Recovered"),
    model.add_state("Deceased")
    ;

    // Setting up parameters
    model.add_param(recovery_rate, "Recovery rate");
    model.add_param(transmission_rate, "Transmission rate"),
    model.add_param(death_rate, "Death rate");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1,2,3);
    virus.set_prob_recovery(&model("Recovery rate"));
    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_prob_death(&model("Death rate"));
    
    model.add_virus(virus, prevalence);

    model.set_name("Susceptible-Infected-Recovered-Deceased (SIRD)");

    return;
   
}

template<typename TSeq>
inline ModelSIRD<TSeq>::ModelSIRD(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    ModelSIRD<TSeq>(
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
