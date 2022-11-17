#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;

[[cpp11::register]]
int init_cpp(SEXP m, int days, int seed) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->init(days, seed);
  
  return 0;
  
}

[[cpp11::register]]
int print_cpp(SEXP m) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->print();
  
  return 0;
  
}

[[cpp11::register]]
int agents_smallworld_cpp(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return 0;
  
}

[[cpp11::register]]
int run_cpp(SEXP m) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->run();
  
  return 0;
  
}

