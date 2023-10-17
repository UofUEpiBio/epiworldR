#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP print_cpp(SEXP m, bool lite) {
  
  external_pointer<Model<>> ptr(m);
  ptr->print(lite);
  
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
  
  external_pointer<Model<>> ptr(m);
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
  
  external_pointer<Model<>> ptr(m);
  ptr->agents_from_edgelist(source, target, size, directed); 
  
  return m;
  
}

[[cpp11::register]]
SEXP run_cpp(SEXP m, int ndays, int seed) {
  
  external_pointer<Model<>> ptr(m);
  ptr->run(ndays, seed);
  
  return m;
  
}

typedef std::function<void(size_t,Model<>*)> funptr;

[[cpp11::register]]
SEXP make_saver_cpp(
  std::string fn,
  bool total_hist,
  bool virus_info,
  bool virus_hist,
  bool tool_info,
  bool tool_hist,
  bool transmission,
  bool transition,
  bool reproductive,
  bool generation
) {
  
  funptr* saver = new funptr(make_save_run<int>(
    fn,
    total_hist,
    virus_info,
    virus_hist,
    tool_info,
    tool_hist,
    transmission,
    transition,
    reproductive, 
    generation
  ));
  
  external_pointer<funptr> sav_ptr(saver);
  
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
  
  external_pointer<Model<>> ptr(m);
  external_pointer<funptr> sav_ptr(saver);
  
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
  
  external_pointer<Model<>> ptr(model);
  ptr->queuing_on();
  return model;
  
}

[[cpp11::register]]
SEXP queuing_off_cpp(
    SEXP model
) {
  
  external_pointer<Model<>> ptr(model);
  ptr->queuing_off();
  return model;
  
}

[[cpp11::register]]
double get_param_cpp(SEXP model, std::string pname) {
  external_pointer<Model<>> ptr(model);
  return static_cast<double>(ptr->get_param(pname));
}

[[cpp11::register]]
SEXP set_param_cpp(SEXP model, std::string pname, double val) {
  
  external_pointer<Model<>> ptr(model);
  ptr->operator()(pname) = val;
  
  return model;
}

[[cpp11::register]]
SEXP set_name_cpp(SEXP model, std::string mname) {
  external_pointer<Model<>> ptr(model);
  ptr->set_name(mname);
  return model;
}

[[cpp11::register]]
std::string get_name_cpp(SEXP model) {
  external_pointer<Model<>> ptr(model);
  return ptr->get_name();
}

[[cpp11::register]]
strings get_states_cpp(
    SEXP model
) {
  
  external_pointer<Model<>> ptr(model);
  return writable::strings(ptr->get_states());
  
}

[[cpp11::register]]
SEXP verbose_on_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  ptr->verbose_on(); 
  return model;
  
}

[[cpp11::register]]
SEXP verbose_off_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  ptr->verbose_off(); 
  return model;
  
}

[[cpp11::register]]
int get_n_viruses_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  return static_cast<int>(ptr->get_n_viruses());
  
}

[[cpp11::register]]
int get_n_tools_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  return static_cast<int>(ptr->get_n_tools());
  
}

[[cpp11::register]]
int get_ndays_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  return static_cast<int>(ptr->get_ndays());
  
}

[[cpp11::register]]
int get_n_replicates_cpp(SEXP model) {
  
  external_pointer<Model<>> ptr(model);
  return static_cast<int>(ptr->get_n_replicates());
  
}

[[cpp11::register]]
int size_cpp(SEXP model) {
  return external_pointer<Model<>>(model)->size();
}

[[cpp11::register]]
SEXP set_agents_data_cpp(SEXP model, SEXP data, int ncols) {

  external_pointer<Model<>> modelptr(model);
  modelptr->set_agents_data(
    REAL(data),
    ncols
  );
  
  return model;
   
}

[[cpp11::register]]
int get_agents_data_ncols_cpp(SEXP model) {
  
  return
    static_cast<int>(
      external_pointer<Model<>>(model)->get_agents_data_ncols()
    );
  
}

[[cpp11::register]]
SEXP get_virus_model_cpp(SEXP model, int virus_pos) {
  external_pointer<Model<>> modelptr(model);
  
  external_pointer<Virus<>> res(
      &modelptr->get_virus(static_cast<size_t>(virus_pos)),
      false
  );
  
  return res;
  
}

[[cpp11::register]]
SEXP get_tool_model_cpp(SEXP model, int tool_pos) {
  external_pointer<Model<>> modelptr(model);
  
  external_pointer<Tool<>> res(
      &modelptr->get_tool(static_cast<size_t>(tool_pos)),
      false
  );
  
  return res;
  
}

[[cpp11::register]]
cpp11::data_frame get_network_cpp(SEXP model) {
    
  external_pointer<Model<>> modelptr(model);
  
  std::vector<int> from;
  std::vector<int> to;

  modelptr->write_edgelist(from, to);

  return cpp11::writable::data_frame({
    "from"_nm = from,
    "to"_nm   = to
  });
    
}

[[cpp11::register]]
SEXP initial_states_cpp(SEXP model, cpp11::doubles proportions) {

  external_pointer<Model<>> modelptr(model);

  std::vector< double > states_vec(proportions.begin(), proportions.end());

  modelptr->initial_states(states_vec, std::vector< int >({}));
  
  return model;
 
}
