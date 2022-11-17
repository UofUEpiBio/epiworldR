
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapSURV(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSURV<>> (a)

[[cpp11::register]]
SEXP ModelSURV_cpp(
  std::string name,
  double prevalence,
  double efficacy_vax,
  double latent_period,         
  double prob_symptoms,         
  double prop_vaccinated,       
  double prop_vax_redux_transm, 
  double infect_period,         
  double prop_vax_redux_infect, 
  double surveillance_prob,     
  double prob_transmission,     
  double prob_death,            
  double prob_noreinfect       
) {
    

  // Creating a pointer to a ModelSIR model
  WrapSURV(ptr)(
    new epiworld::epimodels::ModelSURV<>(
      name,
      prevalence,
      efficacy_vax,
      latent_period,  
      infect_period,  
      prob_symptoms,         
      prop_vaccinated,       
      prop_vax_redux_transm, 
      prop_vax_redux_infect, 
      surveillance_prob,     
      prob_transmission,     
      prob_death,            
      prob_noreinfect   
    )
  );

  return ptr;
}

#undef WrapSURV
