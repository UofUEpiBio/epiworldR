
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;

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
