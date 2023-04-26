
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapTool(a) cpp11::external_pointer<epiworld::Tool<>> (a)

[[cpp11::register]]
SEXP tool_cpp(
    std::string name,
    double susceptibility_reduction,
    double transmission_reduction,
    double recovery_enhancer,
    double death_reduction
    ) {
  
  WrapTool(tool)(new epiworld::Tool<int>(name));
  
  if (susceptibility_reduction > 0)
    tool->set_susceptibility_reduction(susceptibility_reduction);
  
  if (transmission_reduction > 0)
    tool->set_transmission_reduction(transmission_reduction);
  
  if (recovery_enhancer > 0)
    tool->set_recovery_enhancer(recovery_enhancer);
  
  if (death_reduction > 0)
    tool->set_death_reduction(death_reduction);
  

  return tool;
  
}
  
[[cpp11::register]]
int add_tool_cpp(SEXP m, SEXP t, double preval) {
  
  cpp11::external_pointer<epiworld::Model<>>(m)->add_tool(
    *cpp11::external_pointer<epiworld::Tool<>>(t),
    preval
  );
  
  return 0;
}

[[cpp11::register]]
int add_tool_n_cpp(SEXP m, SEXP t, size_t preval) {
  
  cpp11::external_pointer<epiworld::Model<>>(m)->add_tool_n(
      *cpp11::external_pointer<epiworld::Tool<>>(t),
      preval
  );
  
  return 0;
}
  
[[cpp11::register]]
SEXP rm_tool_cpp(SEXP m, size_t tool_pos) {
  cpp11::external_pointer<epiworld::Model<>>(m)->rm_tool(tool_pos);
  return m;
}

// Tool function ---------------------------------------------------------------
[[cpp11::register]]
SEXP tool_fun_logit_cpp(
    integers vars,
    doubles coefs,
    SEXP model
) {
  
  external_pointer<Model<>> mptr(model);  
  
  external_pointer<ToolFun<>> res(
      new ToolFun<>(
          tool_fun_logit(
            as_cpp<std::vector<int>>(vars),
            as_cpp<std::vector<double>>(coefs),
            &(*mptr)
          )
      )
  );
  
  return res;
  
}

// Probability of transmission -------------------------------------------------
[[cpp11::register]]
SEXP set_transmission_reduction_cpp(SEXP tool, double prob) {
  
  WrapTool(toolptr)(tool);
  toolptr->set_transmission_reduction(prob);
  return tool;
  
}

[[cpp11::register]]
SEXP set_transmission_reduction_ptr_cpp(SEXP tool, SEXP model, std::string param) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  
  toolptr->set_transmission_reduction(
      &(mptr->operator()(param))
  );
  
  return tool;
  
}

[[cpp11::register]]
SEXP set_transmission_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  external_pointer<ToolFun<>> tfunptr(tfun);
  
  toolptr->set_transmission_reduction_fun(*tfunptr);
  
  return tool;
  
}

// Probability of recovery -----------------------------------------------------
[[cpp11::register]]
SEXP set_susceptibility_reduction_cpp(SEXP tool, double prob) {
  
  WrapTool(toolptr)(tool);
  toolptr->set_susceptibility_reduction(prob);
  return tool;
  
}

[[cpp11::register]]
SEXP set_susceptibility_reduction_ptr_cpp(SEXP tool, SEXP model, std::string param) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  
  toolptr->set_susceptibility_reduction(
      &(mptr->operator()(param))
  );
  
  return tool;
  
}

[[cpp11::register]]
SEXP set_susceptibility_reduction_cpp(SEXP tool, SEXP model, SEXP tfun) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  external_pointer<ToolFun<>> tfunptr(tfun);
  
  toolptr->set_susceptibility_reduction_fun(*tfunptr);
  
  return tool;
  
}

// Recovery enhancer -----------------------------------------------------------
[[cpp11::register]]
SEXP set_recovery_enhancer_cpp(SEXP tool, double prob) {
  
  WrapTool(toolptr)(tool);
  toolptr->set_recovery_enhancer(prob);
  return tool;
  
}

[[cpp11::register]]
SEXP set_recovery_enhancer_ptr_cpp(SEXP tool, SEXP model, std::string param) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  
  toolptr->set_recovery_enhancer(
      &(mptr->operator()(param))
  );
  
  return tool;
  
}

[[cpp11::register]]
SEXP set_recovery_enhancer_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  external_pointer<ToolFun<>> tfunptr(tfun);
  
  toolptr->set_recovery_enhancer_fun(*tfunptr);
  
  return tool;
  
}

// Death reduction -------------------------------------------------------------
[[cpp11::register]]
SEXP set_death_reduction_cpp(SEXP tool, double prob) {
  
  WrapTool(toolptr)(tool);
  toolptr->set_death_reduction(prob);
  return tool;
  
}

[[cpp11::register]]
SEXP set_death_reduction_ptr_cpp(SEXP tool, SEXP model, std::string param) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  
  toolptr->set_death_reduction(
      &(mptr->operator()(param))
  );
  
  return tool;
  
}

[[cpp11::register]]
SEXP set_death_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  external_pointer<ToolFun<>> tfunptr(tfun);
  
  toolptr->set_death_reduction_fun(*tfunptr);
  
  return tool;
  
}


[[cpp11::register]]
std::string get_name_tool_cpp(SEXP tool) {
  return external_pointer<Tool<>>(tool)->get_name();
}

[[cpp11::register]]
SEXP set_name_tool_cpp(SEXP tool, std::string name) {
  external_pointer<Tool<>>(tool)->set_name(name);
  return tool;
}
  
#undef WrapTool
