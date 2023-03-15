#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;

[[cpp11::register]]
SEXP print_cpp(SEXP m) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->print();
  
  return m;
  
}

[[cpp11::register]]
SEXP agents_smallworld_cpp(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->agents_smallworld(n, k, d, p);
  
  return m;
  
}

[[cpp11::register]]
SEXP run_cpp(SEXP m, int ndays, int seed) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->run(ndays, seed);
  
  return m;
  
}

[[cpp11::register]]
SEXP queuing_on_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->queuing_on();
  return model;
  
}

[[cpp11::register]]
SEXP queuing_off_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->queuing_off();
  return model;
  
}

[[cpp11::register]]
double get_param_cpp(SEXP model, std::string pname) {
  cpp11::external_pointer<Model<>> ptr(model);
  return static_cast<double>(ptr->get_param(pname));
}

[[cpp11::register]]
SEXP set_param_cpp(SEXP model, std::string pname, double val) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->operator()(pname) = val;
  
  return model;
}

[[cpp11::register]]
SEXP set_name_cpp(SEXP model, std::string mname) {
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->set_name(mname);
  return model;
}

[[cpp11::register]]
std::string get_name_cpp(SEXP model) {
  cpp11::external_pointer<Model<>> ptr(model);
  return ptr->get_name();
}

[[cpp11::register]]
cpp11::strings get_state_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::strings(ptr->get_state());
  
}


