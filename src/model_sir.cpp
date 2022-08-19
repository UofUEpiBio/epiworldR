#define printf_epiworld fflush(stdout);Rprintf

#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld.hpp"

using namespace cpp11;

#define WrapSIR(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIR<>> (a)

[[cpp11::register]]
SEXP ModelSIR(
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
int init(SEXP m, int days, int seed) {
  
  WrapSIR(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
};

[[cpp11::register]]
int print(SEXP m) {
  
  WrapSIR(ptr)(m);
  ptr->print();
  
  return 0;
  
}
  
[[cpp11::register]]
int agents_smallworld(
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
int run(SEXP m) {
  
  WrapSIR(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSIR
