#ifndef EPIWORLD_MODELS_SIS_HPP 
#define EPIWORLD_MODELS_SIS_HPP

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Infected-Susceptible (SIS) model
 * 
 * ![Model Diagram](../assets/img/sis.png)
 * 
 * @ingroup basic_compartmental
 * 
 * @param vname std::string Name of the virus
 * @param n Number of agents.
 * @param prevalence epiworld_double Initial share of infected agents
 * @param transmission_prob epiworld_double Probability of transmission per contact
 * @param recovery_rate epiworld_double Daily probability of recovery
 */
class ModelSIS : public Model
{

public:

    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int INFECTED    = 1;

    ModelSIS() = delete;

    ModelSIS(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_prob,
        epiworld_double recovery_rate
    );


    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};

inline ModelSIS::ModelSIS(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double recovery_rate
    )
{

    this->set_population_size(n);

    this->set_name("Susceptible-Infected-Susceptible (SIS)");

    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(recovery_rate, "Recovery rate");

    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Infected", stage::infectious);
    virus.infection(SUSCEPTIBLE, INFECTED, "Transmission probability");
    virus.rate(INFECTED, SUSCEPTIBLE, "Recovery rate").effect(Effect::recovery);

    this->add_virus(virus, assign::share(prevalence));

}

inline ModelInfo ModelSIS::model_info()
{

    using uint_t = epiworld_fast_uint;

    return {"SIS", "Susceptible-Infected-Susceptible",
        "Infected agents become susceptible again when they recover.",
        "Disease",
        {
            ParamInfo::agents(10000),
            ParamInfo::prevalence(0.01),
            ParamInfo::probability("transmission_prob", "Transmission probability",
                0.05, "Probability that a contact with an infected agent infects."),
            ParamInfo::probability("recovery_rate", "Recovery rate", 0.14,
                "Daily probability that an infected agent recovers.")
        },
        [](const ModelArgs & a) -> std::unique_ptr< Model > {
            return std::make_unique< ModelSIS >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("recovery_rate")
            );
        }};
}

#endif
