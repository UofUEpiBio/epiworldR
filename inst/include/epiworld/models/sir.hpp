#ifndef EPIWORLD_SIR_H 
#define EPIWORLD_SIR_H

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIR : public epiworld::Model<TSeq>
{
public:

    ModelSIR() {};

    ModelSIR(
        ModelSIR<TSeq> & model,
        const std::string & vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate
    );

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
    );
    
};

template<typename TSeq>
inline ModelSIR<TSeq>::ModelSIR(
    ModelSIR<TSeq> & model,
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate
    )
{

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Infected", epiworld::default_update_exposed<TSeq>);
    model.add_state("Recovered");

    // Setting up parameters
    model.add_param(recovery_rate, "Recovery rate");
    model.add_param(transmission_rate, "Transmission rate");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1,2,2);
    
    virus.set_prob_recovery(&model("Recovery rate"));
    virus.set_prob_infecting(&model("Transmission rate"));
    
    model.add_virus(virus);

    model.set_name("Susceptible-Infected-Recovered (SIR)");

    return;
   
}

template<typename TSeq>
inline ModelSIR<TSeq>::ModelSIR(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate
    )
{

    ModelSIR<TSeq>(
        *this,
        vname,
        prevalence,
        transmission_rate,
        recovery_rate
        );

    return;

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
