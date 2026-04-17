#ifndef EPIWORLD_MODELS_SEIRD_HPP
#define EPIWORLD_MODELS_SEIRD_HPP

#include "../model-bones.hpp"

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed-Deceased (SEIRD) model
 * 
 * ![Model Diagram](../assets/img/seird.png)
 * 
 * @ingroup death_compartmental
*/
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRD : public Model<TSeq>
{
  
public:
  static const int SUSCEPTIBLE = 0;
  static const int EXPOSED     = 1;
  static const int INFECTED    = 2;
  static const int REMOVED     = 3;
  static const int DECEASED    = 4;
  
  ModelSEIRD() = delete;

  /**
   * @brief Constructor for the SEIRD model.
   * 
   * @param vname Name of the model.
   * @param prevalence Initial prevalence of the disease.
   * @param transmission_rate Transmission rate of the disease.
   * @param avg_incubation_days Average incubation period of the disease.
   * @param recovery_rate Recovery rate of the disease.
   * @param death_rate Death rate of the disease.
   */
  ModelSEIRD(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    epiworld_double death_rate
  );
  
  UpdateFun<TSeq> update_exposed_seir = [](
    Agent<TSeq> * p,
    Model<TSeq> * m
  ) -> void {

    // Getting the virus
    auto v = p->get_virus();

    // Does the agent become infected?
    if (m->runif() < 1.0/(v->get_incubation(m)))
      p->change_state(*m, ModelSEIRD<TSeq>::INFECTED);

    return;
  };

  
  UpdateFun<TSeq> update_infected = [](
    Agent<TSeq> * p, Model<TSeq> * m
  ) -> void {
          
    // Odd: Die, Even: Recover
    epiworld_fast_uint n_events = 0u;

    auto & v = p->get_virus();
      
    // Die
    m->array_double_tmp[n_events++] = 
      v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, *m));
    
    // Recover
    m->array_double_tmp[n_events++] = 
      1.0 - (1.0 - v->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(v, *m));
    
    
#ifdef EPI_DEBUG
    if (n_events == 0u)
    {
      printf_epiworld(
        "[epi-debug] agent %i has 0 possible events!!\n",
        static_cast<int>(p->get_id())
      );
      throw std::logic_error("Zero events in exposed.");
    }
#else
    if (n_events == 0u)
      return;
#endif
    
    
    // Running the roulette
    int which = roulette(n_events, m);
    
    if (which < 0)
      return;
    
    // Which roulette happen?
    if ((which % 2) == 0) // If odd
    {
      
      p->rm_virus(*m, ModelSEIRD<TSeq>::DECEASED);
      
    } else {
      
      p->rm_virus(*m);
      
    }
    
    return ;
      
    
    
    return;
    
  };

  ModelSEIRD<TSeq> & initial_states(
    std::vector< double > proportions_,
    std::vector< int > queue_ = {}
  ) override;

};



template<typename TSeq>
inline ModelSEIRD<TSeq>::ModelSEIRD(
    const std::string & vname,
    epiworld_double prevalence,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    epiworld_double death_rate
)
{
  
  // Adding statuses
  this->add_state("Susceptible", default_update_susceptible<TSeq>);
  this->add_state("Exposed",  this->update_exposed_seir);
  this->add_state("Infected", this->update_infected);
  this->add_state("Removed");
  this->add_state("Deceased");
  
  // Setting up parameters
  this->add_param(transmission_rate, "Transmission rate");
  this->add_param(avg_incubation_days, "Incubation days");
  this->add_param(recovery_rate, "Recovery rate");
  this->add_param(death_rate, "Death rate");
  
  // Preparing the virus -------------------------------------------
  Virus<TSeq> virus(vname, prevalence, true);
  virus.set_state(ModelSEIRD<TSeq>::EXPOSED, ModelSEIRD<TSeq>::REMOVED, ModelSEIRD<TSeq>::DECEASED);
  
  virus.set_prob_infecting("Transmission rate");
  virus.set_incubation("Incubation days");
  virus.set_prob_death("Death rate");
  virus.set_prob_recovery("Recovery rate");
  
  // Adding the tool and the virus
  this->add_virus(virus);
  
  this->set_name("Susceptible-Exposed-Infected-Removed-Deceased (SEIRD)");
  
}

template<typename TSeq>
inline ModelSEIRD<TSeq> & ModelSEIRD<TSeq>::initial_states(
  std::vector< double > proportions_,
  std::vector< int > /**/
) {

  Model<TSeq>::initial_states_fun =
    create_init_function_seird<TSeq>(proportions_)
    ;

  return *this;

}


#endif