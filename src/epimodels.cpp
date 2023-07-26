
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
  double transmission_rate,     
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
      transmission_rate,     
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
        double transmission_rate,
        double incubation_days,
        double recovery_rate
    
    ) {
      
      // Creating a pointer to a ModelSIR model
      WrapSEIR(ptr)(
          new epiworld::epimodels::ModelSEIR<>(
              name,
              prevalence,
              transmission_rate,
              incubation_days,
              recovery_rate
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
    double transmission_rate,
    double recovery_rate
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIS(ptr)(
      new epiworld::epimodels::ModelSIS<>(
          name,
          prevalence,
          transmission_rate,
          recovery_rate
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
    double transmission_rate, 
    double recovery_rate
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIRCONN(ptr)(
      new epiworld::epimodels::ModelSIRCONN<>(
          name,
          n,
          prevalence,
          contact_rate,
          transmission_rate,
          recovery_rate
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
    double transmission_rate,
    double recovery_rate
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIR(ptr)(
      new epiworld::epimodels::ModelSIR<>(
          name,
          prevalence,
          transmission_rate,
          recovery_rate
      )
  );
  
  return ptr;
}
  
#undef WrapSIR

#define WrapSIRD(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIRD<>> (a)

[[cpp11::register]]
SEXP ModelSIRD_cpp(
    std::string name,
    double prevalence,
    double transmission_rate,
    double recovery_rate,
    double death_rate
) {
  
  // Creating a pointer to a ModelSIRD model
  WrapSIRD(ptr)(
      new epiworld::epimodels::ModelSIRD<>(
          name,
          prevalence,
          transmission_rate,
          recovery_rate,
          death_rate
      )
  );
  
  return ptr;
}
  
#undef WrapSIRD
  
#define WrapSISD(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSISD<>> (a)

[[cpp11::register]]
SEXP ModelSISD_cpp(
    std::string name,
    double prevalence,
    double transmission_rate,
    double recovery_rate,
    double death_rate
) {
  
  // Creating a pointer to a ModelSISD model
  WrapSISD(ptr)(
      new epiworld::epimodels::ModelSISD<>(
          name,
          prevalence,
          transmission_rate,
          recovery_rate,
          death_rate
      )
  );
  
  return ptr;
}
  
#undef WrapSISD

#define WrapSURVD(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSURVD<>> (a)

[[cpp11::register]]
SEXP ModelSURVD_cpp(
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
  double transmission_rate,     
  double prob_death,            
  double prob_noreinfect       
) {
    

  // Creating a pointer to a ModelSURVD model
  WrapSURVD(ptr)(
    new epiworld::epimodels::ModelSURVD<>(
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
      transmission_rate,     
      prob_death,            
      prob_noreinfect   
    )
  );

  return ptr;
}

#undef WrapSURVD
  
#define WrapSEIRCONN(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSEIRCONN<>> (a)

[[cpp11::register]]
SEXP ModelSEIRCONN_cpp(
    std::string name,
    unsigned int n,
    double prevalence,
    double contact_rate,
    double transmission_rate,
    double incubation_days,
    double recovery_rate
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSEIRCONN(ptr)(
      new epiworld::epimodels::ModelSEIRCONN<>(
          name,
          n,
          prevalence,
          contact_rate,
          transmission_rate,
          incubation_days,
          recovery_rate
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
  double recovery_rate,
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
        recovery_rate,
        prevalence
    )
  );
  
  return ptr;
  
}


[[cpp11::register]]
SEXP ModelDiffNet_cpp(
  std::string name,
  double prevalence,
  double prob_adopt,
  bool normalize_exposure,
  SEXP data,
  int data_ncols,
  std::vector< int > data_cols,
  std::vector<double> params
) {

  // Maps data_cols to size_t
  std::vector< size_t > data_cols_s;
  for (auto i : data_cols)
    data_cols_s.push_back(static_cast<size_t>(i));
  
  cpp11::external_pointer<epiworld::epimodels::ModelDiffNet<>> ptr(
    new epiworld::epimodels::ModelDiffNet<>(
      name,
      prevalence,
      prob_adopt,
      normalize_exposure,
      &(REAL(data)[0u]),
      data_ncols,
      data_cols_s,
      params
    )
  );
  
  return ptr;
  
}

