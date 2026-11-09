#ifndef EPIWORLD_MODELS_SISD_HPP 
#define EPIWORLD_MODELS_SISD_HPP

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Infected-Susceptible-Deceased (SISD) model
 * 
 * ![Model Diagram](../assets/img/sisd.png)
 * 
 * @ingroup death_compartmental
 * 
 * @param vname std::string Name of the virus
 * @param n Number of agents.
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery_rate rate of the immune system
 * @param inital_death epiworld_double Initial death_rate of the immune system
 */
class ModelSISD : public Model
{

public:

    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int INFECTED    = 1;
    static constexpr int DECEASED    = 2;

    ModelSISD() = delete;

    ModelSISD(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_prob,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );


    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};

inline ModelSISD::ModelSISD(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    this->set_population_size(n);

    this->set_name("Susceptible-Infected-Susceptible-Deceased (SISD)");

    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(death_rate, "Death rate");

    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Infected", stage::infectious);
    virus.add_stage("Deceased");
    virus.set_contact(DECEASED, 0.0);
    virus.infection(SUSCEPTIBLE, INFECTED, "Transmission probability");
    virus.rate(INFECTED, SUSCEPTIBLE, "Recovery rate").effect(Effect::recovery);
    virus.rate(INFECTED, DECEASED, "Death rate").effect(Effect::death);

    this->add_virus(virus, assign::share(prevalence));

}

inline ModelInfo ModelSISD::model_info()
{

    using uint_t = epiworld_fast_uint;

    return {"SISD", "Susceptible-Infected-Susceptible-Deceased",
        "SIS, where infected agents may also die.",
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
            return std::make_unique< ModelSISD >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("recovery_rate"), a("death_rate")
            );
        }};
}

#endif
