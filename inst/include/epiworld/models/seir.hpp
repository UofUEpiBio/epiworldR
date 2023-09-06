#ifndef EPIWORLD_MODELS_SEIR_HPP
#define EPIWORLD_MODELS_SEIR_HPP

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence epiworld_double Initial prevalence the immune system
 * @param transmission_rate epiworld_double Transmission rate of the virus
 * @param avg_incubation_days epiworld_double Average incubation days of the virus.
 * @param recovery_rate epiworld_double Recovery rate of the virus.
 */
template<typename TSeq = int>
class ModelSEIR : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int REMOVED     = 3;

public:

    ModelSEIR() {};

    ModelSEIR(
        ModelSEIR<TSeq> & model,
        std::string vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate
    );

    ModelSEIR(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate
    );
    
    epiworld::UpdateFun<TSeq> update_exposed_seir = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Model<TSeq> * m
    ) -> void {

        // Getting the virus
        auto v = p->get_virus(0);

        // Does the agent become infected?
        if (m->runif() < 1.0/(v->get_incubation(m)))
            p->change_state(m, ModelSEIR<TSeq>::INFECTED);

        return;    
    };
      

    epiworld::UpdateFun<TSeq> update_infected_seir = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Model<TSeq> * m
    ) -> void {
        // Does the agent recover?
        if (m->runif() < (m->par("Recovery rate")))
            p->rm_virus(0, m);

        return;    
    };

};


template<typename TSeq>
inline ModelSEIR<TSeq>::ModelSEIR(
    ModelSEIR<TSeq> & model,
    std::string vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate
    )
{

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Exposed", model.update_exposed_seir);
    model.add_state("Infected", model.update_infected_seir);
    model.add_state("Removed");

    // Setting up parameters
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(avg_incubation_days, "Incubation days");
    model.add_param(recovery_rate, "Recovery rate");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(ModelSEIR<TSeq>::EXPOSED, ModelSEIR<TSeq>::REMOVED, ModelSEIR<TSeq>::REMOVED);

    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_incubation(&model("Incubation days"));
    
    // Adding the tool and the virus
    model.add_virus(virus, prevalence);
    
    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR)");

    return;
   
}

template<typename TSeq>
inline ModelSEIR<TSeq>::ModelSEIR(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate
    )
{

    ModelSEIR<TSeq>(
        *this,
        vname,
        prevalence,
        transmission_rate,
        avg_incubation_days,
        recovery_rate
        );

    return;

}



#endif
