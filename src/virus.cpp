
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;


// Model definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/models

#define WrapVirus(a) external_pointer<Virus<>> (a)

[[cpp11::register]]
SEXP virus_cpp(
    std::string name,
    double prob_infecting,
    double prob_recovery,
    double prob_death,
    double post_immunity
    ) {
  
  WrapVirus(virus)(new epiworld::Virus<int>(name));
  
  virus->set_prob_infecting(prob_infecting);
  virus->set_prob_recovery(prob_recovery);
  virus->set_prob_death(prob_death);
  
  if (post_immunity > 0.0)
    virus->set_post_immunity(post_immunity);
  
  return virus;
  
}

[[cpp11::register]]
SEXP virus_set_state_cpp(
  SEXP v,
  size_t init,
  size_t end,
  size_t removed
) {
  
  WrapVirus(vptr)(v);
  vptr->set_state(
    init, end, removed
  );
  
  return v;
  
}

[[cpp11::register]]
SEXP add_virus_cpp(SEXP m, SEXP v, double preval) {
  
  external_pointer<epiworld::Model<>>(m)->add_virus(
    *external_pointer<epiworld::Virus<>>(v),
    preval
  );
  
  return m;
}

[[cpp11::register]]
SEXP add_virus_n_cpp(SEXP m, SEXP v, size_t preval) {
  
  external_pointer<Model<>>(m)->add_virus_n(
      *external_pointer<Virus<>>(v),
      preval
  );
  
  return m;
}

[[cpp11::register]]
SEXP rm_virus_cpp(SEXP m, size_t virus_pos) {
  
  external_pointer<epiworld::Model<>>(m)->rm_virus(virus_pos);
  return m;
  
}

[[cpp11::register]]
SEXP print_virus_cpp(SEXP v) {
  
  WrapVirus(vptr)(v);
  vptr->print();
  return v;
  
}

// Virus function --------------------------------------------------------------
[[cpp11::register]]
SEXP virus_fun_logit_cpp(
  integers vars,
  doubles coefs,
  SEXP model
  ) {
  
  external_pointer<Model<>> mptr(model);  
  
  external_pointer<VirusFun<>> res(
      new VirusFun<>(
          virus_fun_logit(
            as_cpp<std::vector<int>>(vars),
            as_cpp<std::vector<double>>(coefs),
            &(*mptr)
          )
      )
  );
  
  return res;
  
}
 
// Probability of infection ----------------------------------------------------
[[cpp11::register]]
SEXP set_prob_infecting_cpp(SEXP virus, double prob) {
  
  WrapVirus(vptr)(virus);
  vptr->set_prob_infecting(prob);
  return virus;
  
}
  
[[cpp11::register]]
SEXP set_prob_infecting_ptr_cpp(SEXP virus, SEXP model, std::string param) {
  
  WrapVirus(vptr)(virus);
  external_pointer<Model<>> mptr(model);
  
  vptr->set_prob_infecting(
    &(mptr->operator()(param))
  );
  
  return virus;
  
}
  
[[cpp11::register]]
SEXP set_prob_infecting_fun_cpp(SEXP virus, SEXP model, SEXP vfun) {
  
  WrapVirus(vptr)(virus);
  external_pointer<Model<>> mptr(model);
  external_pointer<VirusFun<>> vfunptr(vfun);
  
  vptr->set_prob_infecting_fun(*vfunptr);
  
  return virus;
  
}
 
// Probability of recovery -----------------------------------------------------
[[cpp11::register]]
SEXP set_prob_recovery_cpp(SEXP virus, double prob) {
 
 WrapVirus(vptr)(virus);
 vptr->set_prob_recovery(prob);
 return virus;
 
}

[[cpp11::register]]
SEXP set_prob_recovery_ptr_cpp(SEXP virus, SEXP model, std::string param) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 
 vptr->set_prob_recovery(
     &(mptr->operator()(param))
 );
 
 return virus;
 
}

[[cpp11::register]]
SEXP set_prob_recovery_fun_cpp(SEXP virus, SEXP model, SEXP vfun) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 external_pointer<VirusFun<>> vfunptr(vfun);
 
 vptr->set_prob_recovery_fun(*vfunptr);
 
 return virus;
 
}
 
// Probability of death --------------------------------------------------------
[[cpp11::register]]
SEXP set_prob_death_cpp(SEXP virus, double prob) {
 
 WrapVirus(vptr)(virus);
 vptr->set_prob_death(prob);
 return virus;
 
}

[[cpp11::register]]
SEXP set_prob_death_ptr_cpp(SEXP virus, SEXP model, std::string param) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 
 vptr->set_prob_death(
     &(mptr->operator()(param))
 );
 
 return virus;
 
}

[[cpp11::register]]
SEXP set_prob_death_fun_cpp(SEXP virus, SEXP model, SEXP vfun) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 external_pointer<VirusFun<>> vfunptr(vfun);
 
 vptr->set_prob_death_fun(*vfunptr);

 return virus;
 
}

// Incubation period ----------------------------------------------------------
[[cpp11::register]]
SEXP set_incubation_cpp(SEXP virus, double prob) {
 
 WrapVirus(vptr)(virus);
 vptr->set_incubation(prob);
 return virus;
 
}

[[cpp11::register]]
SEXP set_incubation_ptr_cpp(SEXP virus, SEXP model, std::string param) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 
 vptr->set_incubation(
     &(mptr->operator()(param))
 );
 
 return virus;
 
}

[[cpp11::register]]
SEXP set_incubation_fun_cpp(SEXP virus, SEXP model, SEXP vfun) {
 
 WrapVirus(vptr)(virus);
 external_pointer<Model<>> mptr(model);
 external_pointer<VirusFun<>> vfunptr(vfun);
 
 vptr->set_incubation_fun(*vfunptr);
 
 return virus;
 
}

[[cpp11::register]]
std::string get_name_virus_cpp(SEXP virus) {
  return external_pointer<Virus<>>(virus)->get_name();
}

[[cpp11::register]]
SEXP set_name_virus_cpp(SEXP virus, std::string name) {
  external_pointer<Virus<>>(virus)->set_name(name);
  return virus;
}
    
#undef WrapVirus
