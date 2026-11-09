#ifndef EPIWORLD_MODELS_SEIR_HPP
#define EPIWORLD_MODELS_SEIR_HPP

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * ![Model Diagram](../assets/img/seir.png)
 * 
 * @ingroup basic_compartmental
 *
 * @param model A Model object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param n Number of agents.
 * @param prevalence epiworld_double Initial prevalence the immune system
 * @param transmission_prob epiworld_double Probability that a contact with an
 * infected agent infects.
 * @param latent_period epiworld_double Mean days in Exposed (carrying the
 * virus, not yet infectious).
 * @param recovery_rate epiworld_double Recovery rate of the virus.
 */
class ModelSEIR : public Model
{

public:
    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int EXPOSED     = 1;
    static constexpr int INFECTED    = 2;
    static constexpr int REMOVED     = 3;

    ModelSEIR() = delete;

    ModelSEIR(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double transmission_prob,
        epiworld_double latent_period,
        epiworld_double recovery_rate
    );



    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};


inline ModelSEIR::ModelSEIR(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double latent_period,
    epiworld_double recovery_rate
    )
{

    this->set_population_size(n);

    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(latent_period, "Latent period");
    this->add_param(recovery_rate, "Recovery rate");

    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Exposed", stage::carrier);
    virus.add_stage("Infected", stage::infectious);
    virus.add_stage("Removed");
    virus.infection(SUSCEPTIBLE, EXPOSED, "Transmission probability");
    virus.after(EXPOSED, INFECTED, Duration::geometric("Latent period"));
    virus.rate(INFECTED, REMOVED, "Recovery rate").effect(Effect::recovery);

    this->add_virus(virus, assign::share(prevalence));

    this->set_name("Susceptible-Exposed-Infected-Removed (SEIR)");

}


inline ModelInfo ModelSEIR::model_info()
{

    using uint_t = epiworld_fast_uint;

    return {"SEIR", "Susceptible-Exposed-Infected-Recovered",
        "Infected agents go through a latent (exposed, not infectious) period "
        "before they become infectious.",
        "Disease",
        {
            ParamInfo::agents(10000),
            ParamInfo::prevalence(0.01),
            ParamInfo::probability("transmission_prob", "Transmission probability",
                0.05, "Probability that a contact with an infected agent infects."),
            ParamInfo::duration("latent_period", "Latent period", 7.0,
                "Mean days in the exposed state (infected, not yet infectious)."),
            ParamInfo::probability("recovery_rate", "Recovery rate", 0.14,
                "Daily probability that an infected agent recovers.")
        },
        [](const ModelArgs & a) -> std::unique_ptr< Model > {
            return std::make_unique< ModelSEIR >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("latent_period"), a("recovery_rate")
            );
        }};
}

#endif
