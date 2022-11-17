
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapSIRCONN(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIRCONN<>> (a)

[[cpp11::register]]
SEXP ModelSIRCONN_cpp(
  std::string name,
  unsigned int n,

  double prevalence,
  double reproductive_number,
  double prob_transmission, 
  double prob_recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIRCONN(ptr)(
    new epiworld::epimodels::ModelSIRCONN<>(
      name,
      n,
      prevalence,
      reproductive_number,
      prob_transmission,
      prob_recovery
    )
  );
  
  return ptr;
}

#undef WrapSIRCONN
