
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/matrix.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/models

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

#define WrapSEIR(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSEIR<>> (a)

[[cpp11::register]]
    SEXP ModelSEIR_cpp(
        std::string name,
        double prevalence,
        double infectiousness,
        double incubation_days,
        double recovery
    
    ) {
      
      // Creating a pointer to a ModelSIR model
      WrapSEIR(ptr)(
          new epiworld::epimodels::ModelSEIR<>(
              name,
              prevalence,
              infectiousness,
              incubation_days,
              recovery
          )
      );
      
      return ptr;
    }
  
#undef WrapSEIR
  
#define WrapSIS(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIS<>> (a)

[[cpp11::register]]
SEXP ModelSIS_cpp(
    std::string name,
    double prevalence,
    double infectiousness,
    double recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIS(ptr)(
      new epiworld::epimodels::ModelSIS<>(
          name,
          prevalence,
          infectiousness,
          recovery
      )
  );
  
  
  return ptr;
}
  
  
#undef WrapSIS
  
  
#define WrapSIRCONN(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIRCONN<>> (a)

[[cpp11::register]]
SEXP ModelSIRCONN_cpp(
    std::string name,
    unsigned int n,
    
    double prevalence,
    double contact_rate,
    double prob_transmission, 
    double prob_recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIRCONN(ptr)(
      new epiworld::epimodels::ModelSIRCONN<>(
          name,
          n,
          prevalence,
          contact_rate,
          prob_transmission,
          prob_recovery
      )
  );
  
  return ptr;
}
  
#undef WrapSIRCONN

  
#define WrapSIR(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIR<>> (a)

[[cpp11::register]]
SEXP ModelSIR_cpp(
    std::string name,
    double prevalence,
    double infectiousness,
    double recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIR(ptr)(
      new epiworld::epimodels::ModelSIR<>(
          name,
          prevalence,
          infectiousness,
          recovery
      )
  );
  
  return ptr;
}
  
#undef WrapSIR

  
#define WrapSEIRCONN(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSEIRCONN<>> (a)

[[cpp11::register]]
SEXP ModelSEIRCONN_cpp(
    std::string name,
    unsigned int n,
    double prevalence,
    double contact_rate,
    double prob_transmission,
    double incubation_days,
    double prob_recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSEIRCONN(ptr)(
      new epiworld::epimodels::ModelSEIRCONN<>(
          name,
          n,
          prevalence,
          contact_rate,
          prob_transmission,
          incubation_days,
          prob_recovery
      )
  );
  
  return ptr;
}
  
  
#undef WrapSEIRCONN
  
[[cpp11::register]]
SEXP ModelSIRLogit_cpp(
  std::string vname,
  SEXP data,
  int ncols,
  std::vector< double > coefs_infect,
  std::vector< double > coefs_recover,
  std::vector< int > coef_infect_cols,
  std::vector< int > coef_recover_cols,
  double prob_infection,
  double prob_recovery,
  double prevalence
) {
  
  std::vector< size_t > cinfect;
  std::vector< size_t > crecover;
  
  for (auto i : coef_infect_cols)
    cinfect.push_back(static_cast<size_t>(i));
  
  for (auto i : coef_recover_cols)
    crecover.push_back(static_cast<size_t>(i));
  
  cpp11::external_pointer<epiworld::epimodels::ModelSIRLogit<>> ptr(
    new epiworld::epimodels::ModelSIRLogit<>(
        vname,
        REAL(data),
        static_cast<size_t>(ncols),
        coefs_infect,
        coefs_recover,
        cinfect,
        crecover,
        prob_infection,
        prob_recovery,
        prevalence
    )
  );
  
  return ptr;
  
}
