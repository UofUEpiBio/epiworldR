
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

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
