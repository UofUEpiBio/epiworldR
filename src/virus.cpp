
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapVirus(a) \
  cpp11::external_pointer<epiworld::Virus<>> (a)

[[cpp11::register]]
SEXP virus_cpp(
    std::string name,
    double post_immunity,
    double prob_infecting,
    double prob_recovery,
    double prob_death
    ) {
  
  WrapVirus(virus)(new epiworld::Virus<int>(name));
  
  if (post_immunity > 0.0)
    virus->set_post_immunity(post_immunity);
  
  virus->set_prob_infecting(prob_infecting);
  virus->set_prob_recovery(prob_recovery);
  virus->set_prob_death(prob_death);

  return virus;
  
}
  
[[cpp11::register]]
int virus_set_status_cpp(
  SEXP v,
  size_t init,
  size_t end,
  size_t removed
) {
  
  WrapVirus(vptr)(v);
  vptr->set_status(
    init, end, removed
  );
  
  return 0;
  
}

[[cpp11::register]]
int add_virus_cpp(SEXP m, SEXP v, double preval) {
  
  cpp11::external_pointer<epiworld::Model<>>(m)->add_virus(
    *cpp11::external_pointer<epiworld::Virus<>>(v),
    preval
  );
  
  return 0;
}

[[cpp11::register]]
int add_virus_n_cpp(SEXP m, SEXP v, size_t preval) {
  
  cpp11::external_pointer<epiworld::Model<>>(m)->add_virus_n(
      *cpp11::external_pointer<epiworld::Virus<>>(v),
      preval
  );
  
  return 0;
}
  
#undef WrapVirus
