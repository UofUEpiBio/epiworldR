#ifndef EPIWORLD_SIRD_H 
#define EPIWORLD_SIRD_H

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Infected-Removed-Deceased (SIRD) model
 * 
 * ![Model Diagram](../assets/img/sird.png)
 * 
 * @ingroup death_compartmental
 */
class ModelSIRD : public Model
{
public:

    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int INFECTED    = 1;
    static constexpr int RECOVERED   = 2;
    static constexpr int DECEASED    = 3;

    ModelSIRD() = delete;

    /**
     * @brief Constructs a new SIRD model with the given parameters.
     * 
     * @param vname The name of the vertex associated with this model.
 * @param n Number of agents.
     * @param prevalence The initial prevalence of the disease in the population.
     * @param transmission_prob Probability that a contact with an infected agent infects.
     * @param recovery_rate The rate at which infected individuals recover and become immune.
     * @param death_rate The rate at which infected individuals die.
     */
    ///@{
    ModelSIRD(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_prob,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );
    ///@}

    


    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};

inline ModelSIRD::ModelSIRD(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double recovery_rate, 
    epiworld_double death_rate
    )
{

    this->set_population_size(n);

    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(death_rate, "Death rate");

    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Infected", stage::infectious);
    virus.add_stage("Recovered");
    virus.add_stage("Deceased");
    virus.set_contact(DECEASED, 0.0);
    virus.infection(SUSCEPTIBLE, INFECTED, "Transmission probability");
    virus.rate(INFECTED, RECOVERED, "Recovery rate").effect(Effect::recovery);
    virus.rate(INFECTED, DECEASED, "Death rate").effect(Effect::death);

    this->add_virus(virus, assign::share(prevalence));

    this->set_name("Susceptible-Infected-Recovered-Deceased (SIRD)");

   
}


inline ModelInfo ModelSIRD::model_info()
{

    using uint_t = epiworld_fast_uint;

    return {"SIRD", "Susceptible-Infected-Recovered-Deceased",
        "SIR, where infected agents may also die.",
        "Disease",
        {
            ParamInfo::agents(10000),
            ParamInfo::prevalence(0.01),
            ParamInfo::probability("transmission_prob", "Transmission probability",
                0.05, "Probability that a contact with an infected agent infects."),
            ParamInfo::probability("recovery_rate", "Recovery rate", 0.14,
                "Daily probability that an infected agent recovers."),
            ParamInfo::probability("death_rate", "Death rate", 0.01,
                "Daily probability that an infected agent dies.")
        },
        [](const ModelArgs & a) -> std::unique_ptr< Model > {
            return std::make_unique< ModelSIRD >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("recovery_rate"), a("death_rate")
            );
        }};
}

#endif
