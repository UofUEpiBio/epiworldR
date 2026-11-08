#ifndef EPIWORLD_MODELS_SEIRD_HPP
#define EPIWORLD_MODELS_SEIRD_HPP

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Exposed-Infected-Removed-Deceased (SEIRD) model
 * 
 * ![Model Diagram](../assets/img/seird.png)
 * 
 * @ingroup death_compartmental
*/
class ModelSEIRD : public Model
{
  
public:
  static constexpr int SUSCEPTIBLE = 0;
  static constexpr int EXPOSED     = 1;
  static constexpr int INFECTED    = 2;
  static constexpr int REMOVED     = 3;
  static constexpr int DECEASED    = 4;
  
  ModelSEIRD() = delete;

  /**
   * @brief Constructor for the SEIRD model.
   * 
   * @param vname Name of the model.
 * @param n Number of agents.
   * @param prevalence Initial prevalence of the disease.
   * @param transmission_prob Probability that a contact with an infected
   * agent infects.
   * @param latent_period Mean days in Exposed (carrying the virus, not yet
   * infectious).
   * @param recovery_rate Recovery rate of the disease.
   * @param death_rate Death rate of the disease.
   */
  ModelSEIRD(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double latent_period,
    epiworld_double recovery_rate,
    epiworld_double death_rate
  );
  


    /// This model's entry in a `ModelRegistry`: arguments and builder
    static ModelInfo model_info();
};



inline ModelSEIRD::ModelSEIRD(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double transmission_prob,
    epiworld_double latent_period,
    epiworld_double recovery_rate,
    epiworld_double death_rate
)
{

    this->set_population_size(n);
  
  this->add_param(transmission_prob, "Transmission probability");
  this->add_param(latent_period, "Latent period");
  this->add_param(recovery_rate, "Recovery rate");
  this->add_param(death_rate, "Death rate");

  Virus virus(vname);
  virus.add_stage("Susceptible");
  virus.add_stage("Exposed", stage::carrier);
  virus.add_stage("Infected", stage::infectious);
  virus.add_stage("Removed");
  virus.add_stage("Deceased");
  virus.set_contact(DECEASED, 0.0);
  virus.infection(SUSCEPTIBLE, EXPOSED, "Transmission probability");
  virus.after(EXPOSED, INFECTED, Duration::geometric("Latent period"));
  virus.rate(INFECTED, REMOVED, "Recovery rate").effect(Effect::recovery);
  virus.rate(INFECTED, DECEASED, "Death rate").effect(Effect::death);

  this->add_virus(virus, assign::share(prevalence));

  this->set_name("Susceptible-Exposed-Infected-Removed-Deceased (SEIRD)");
  
}



inline ModelInfo ModelSEIRD::model_info()
{

    using uint_t = epiworld_fast_uint;

    return {"SEIRD", "Susceptible-Exposed-Infected-Recovered-Deceased",
        "SEIR, where infected agents may also die.",
        "Disease",
        {
            ParamInfo::agents(10000),
            ParamInfo::prevalence(0.01),
            ParamInfo::probability("transmission_prob", "Transmission probability",
                0.05, "Probability that a contact with an infected agent infects."),
            ParamInfo::duration("latent_period", "Latent period", 7.0,
                "Mean days in the exposed state (infected, not yet infectious)."),
            ParamInfo::probability("recovery_rate", "Recovery rate", 0.14,
                "Daily probability that an infected agent recovers."),
            ParamInfo::probability("death_rate", "Death rate", 0.01,
                "Daily probability that an infected agent dies.")
        },
        [](const ModelArgs & a) -> std::unique_ptr< Model > {
            return std::make_unique< ModelSEIRD >(
                a.virus_name(), static_cast< uint_t >(a("n")), a("prevalence"),
                a("transmission_prob"), a("latent_period"), a("recovery_rate"),
                a("death_rate")
            );
        }};
}

#endif