
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

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
