
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.hpp"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapSIRCONN(a) \
  cpp11::external_pointer<epiworld::epimodels::ModelSIRCONN<>> (a)

[[cpp11::register]]
<<<<<<< HEAD
SEXP ModelSIRCONN(
  std::string name,
=======
SEXP ModelSIRCONN_cpp(
  std::string name,
  unsigned int n,
>>>>>>> c4b483bb37e81678bcb5a10721bf8633e1d05332
  double prevalence,
  double reproductive_number,
  double prob_transmission, 
  double prob_recovery
) {
  
  // Creating a pointer to a ModelSIR model
  WrapSIRCONN(ptr)(
    new epiworld::epimodels::ModelSIRCONN<>(
      name,
<<<<<<< HEAD
=======
      n,
>>>>>>> c4b483bb37e81678bcb5a10721bf8633e1d05332
      prevalence,
      reproductive_number,
      prob_transmission,
      prob_recovery
    )
  );
  
  return ptr;
}

[[cpp11::register]]
int init_sirconn(SEXP m, int days, int seed) {
  
  WrapSIRCONN(ptr)(m);
  ptr->init(days, seed);
  
  return 0;
  
}
  
[[cpp11::register]]
int print_sirconn(SEXP m) {
  
  WrapSIRCONN(ptr)(m);
  ptr->print();
  
  return 0;
  
}
<<<<<<< HEAD
  
[[cpp11::register]]
int agents_smallworld_sirconn(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

  ) {
  
  WrapSIRCONN(ptr)(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}
=======
>>>>>>> c4b483bb37e81678bcb5a10721bf8633e1d05332

[[cpp11::register]]
int run_sirconn(SEXP m) {
  
  WrapSIRCONN(ptr)(m);
  ptr->run();
  
  return 0;
  
}

#undef WrapSIRCONN
