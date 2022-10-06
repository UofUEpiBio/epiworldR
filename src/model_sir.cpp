
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

[[cpp11::register]]
int init_sir(SEXP m, int days, int seed) {
  
  WrapSIR(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
}

[[cpp11::register]]
int print_sir(SEXP m) {
  
  WrapSIR(ptr)(m);
  ptr->print();
  
  return 0;
  
}
  
[[cpp11::register]]
int agents_smallworld_sir(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

  ) {
  
  WrapSIR(ptr)(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}

[[cpp11::register]]
int run_sir(SEXP m) {
  
  WrapSIR(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSIR
