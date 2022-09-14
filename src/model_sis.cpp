
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.hpp"

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

[[cpp11::register]]
int init_sis(SEXP m, int days, int seed) {
  
  WrapSIS(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
}
  
[[cpp11::register]]
int print_sis(SEXP m) {
  
  WrapSIS(ptr)(m);
  ptr->print();
  
  return 0;
  
}
  
[[cpp11::register]]
int agents_smallworld_sis(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

  ) {
  
  WrapSIS(ptr)(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}

[[cpp11::register]]
int run_sis(SEXP m) {
  
  WrapSIS(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSIS
