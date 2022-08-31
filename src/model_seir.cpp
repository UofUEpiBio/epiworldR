
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.hpp"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapSEIR(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSEIR<>> (a)

[[cpp11::register]]
SEXP ModelSEIR(
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

[[cpp11::register]]
int init_seir(SEXP m, int days, int seed) {
  
  WrapSEIR(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
}
  
[[cpp11::register]]
int print_seir(SEXP m) {
  
  WrapSEIR(ptr)(m);
  ptr->print();
  
  return 0;
  
}
  
[[cpp11::register]]
int agents_smallworld_seir(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

  ) {
  
  WrapSEIR(ptr)(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}

[[cpp11::register]]
int run_seir(SEXP m) {
  
  WrapSEIR(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSEIR
