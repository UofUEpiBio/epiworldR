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

[[cpp11::register]]
int queuing_on_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->queuing_on();
  return 0;
  
}

[[cpp11::register]]
int queuing_off_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->queuing_off();
  return 0;
  
}

[[cpp11::register]]
double get_param_cpp(SEXP model, std::string pname) {
  cpp11::external_pointer<Model<>> ptr(model);
  return static_cast<double>(ptr->get_param(pname));
}

[[cpp11::register]]
int set_param_cpp(SEXP model, std::string pname, double val) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->operator()(pname) = val;
  
  return 0;
}

[[cpp11::register]]
int set_name_cpp(SEXP model, std::string mname) {
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->set_name(mname);
  return 0;
}

[[cpp11::register]]
std::string get_name_cpp(SEXP model) {
  cpp11::external_pointer<Model<>> ptr(model);
  return ptr->get_name();
}

[[cpp11::register]]
cpp11::strings get_status_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::strings(ptr->get_status());
  
}


