#ifndef EPIWORLD_SIRD_H 
#define EPIWORLD_SIRD_H

/**
 * @brief Template for a Susceptible-Infected-Removed-Deceased (SIRD) model
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRD : public epiworld::Model<TSeq>
{
public:

    ModelSIRD() {};

    
    /**
     * @brief Constructs a new SIRD model with the given parameters.
     * 
     * @param model The SIRD model to copy from.
     * @param vname The name of the vertex associated with this model.
     * @param prevalence The initial prevalence of the disease in the population.
     * @param transmission_rate The rate at which the disease spreads from infected to susceptible individuals.
     * @param recovery_rate The rate at which infected individuals recover and become immune.
     * @param death_rate The rate at which infected individuals die.
     */
    ///@{
    ModelSIRD(
        ModelSIRD<TSeq> & model,
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );

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
    ModelSIRD<TSeq> & model,
    const std::string & vname,
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
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,2,3);
    virus.set_prob_recovery(&model("Recovery rate"));
    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_prob_death(&model("Death rate"));
    
    model.add_virus(virus);

    model.set_name("Susceptible-Infected-Recovered-Deceased (SIRD)");

    return;
   
}

template<typename TSeq>
inline ModelSIRD<TSeq>::ModelSIRD(
    const std::string & vname,
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
