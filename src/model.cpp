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
SEXP agents_from_edgelist_cpp(
  SEXP m,
  const std::vector<int> & source,
  const std::vector<int> & target,
  int size,
  bool directed
) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->agents_from_edgelist(source, target, size, directed); 
  
  return m;
  
}

[[cpp11::register]]
SEXP run_cpp(SEXP m, int ndays, int seed) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  ptr->run(ndays, seed);
  
  return m;
  
}

typedef std::function<void(size_t,Model<>*)> funptr;

[[cpp11::register]]
SEXP make_saver_cpp(
  std::string fn,
  bool total_hist,
  bool variant_info,
  bool variant_hist,
  bool tool_info,
  bool tool_hist,
  bool transmission,
  bool transition,
  bool reproductive,
  bool generation
) {
  
  funptr* saver = new funptr(make_save_run<>(
    fn,
    total_hist,
    variant_info,
    variant_hist,
    tool_info,
    tool_hist,
    transmission,
    transition,
    reproductive, 
    generation
  ));
  
  cpp11::external_pointer<funptr> sav_ptr(saver);
  
  return sav_ptr;
  
}

[[cpp11::register]]
SEXP run_multiple_cpp(
    SEXP m,
    int ndays,
    int nsims,
    int seed,
    SEXP saver,
    bool reset,
    bool verbose,
    int nthreads
) {
  
  cpp11::external_pointer<Model<>> ptr(m);
  cpp11::external_pointer<funptr> sav_ptr(saver);
  
  ptr->run_multiple(
    static_cast< epiworld_fast_uint >(ndays),
    static_cast< epiworld_fast_uint >(nsims),
    seed,
    *sav_ptr,
    reset,
    verbose,
    nthreads
  );
  
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

[[cpp11::register]]
SEXP verbose_on_cpp(SEXP model) {
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->verbose_on(); 
  return model;
}

[[cpp11::register]]
SEXP verbose_off_cpp(SEXP model) {
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->verbose_off(); 
  return model;
}

