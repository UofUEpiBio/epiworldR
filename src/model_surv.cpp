
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.hpp"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapSURV(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSURV<>> (a)

[[cpp11::register]]
SEXP ModelSURV(
  std::string name,
  double prevalence,
  double infectiousness,
  double recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSURV(ptr)(
    new epiworld::epimodels::ModelSURV<>(
      name,
      prevalence,
      infectiousness,
      recovery
    )
  );
  
  
  return ptr;
}

[[cpp11::register]]
int init_surv(SEXP m, int days, int seed) {
  
  WrapSURV(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
}
  
[[cpp11::register]]
int print_surv(SEXP m) {
  
  WrapSURV(ptr)(m);
  ptr->print();
  
  return 0;
  
}
  
[[cpp11::register]]
int agents_smallworld_surv(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

  ) {
  
  WrapSURV(ptr)(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}

[[cpp11::register]]
int run_surv(SEXP m) {
  
  WrapSURV(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSURV
