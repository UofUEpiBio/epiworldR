
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;


// Model definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/models

#define WrapTool(a) cpp11::external_pointer<epiworld::Tool<>> (a)

[[cpp11::register]]
SEXP tool_cpp(
    std::string name,
    double prevalence,
    bool as_proportion,
    double susceptibility_reduction,
    double transmission_reduction,
    double recovery_enhancer,
    double death_reduction
    ) {
  
  WrapTool(tool)(new epiworld::Tool<int>(
    name,
    prevalence,
    as_proportion
    ));
  
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
int add_tool_cpp(SEXP m, SEXP t) {
  
  cpp11::external_pointer<epiworld::Model<>>(m)->add_tool(
    *cpp11::external_pointer<epiworld::Tool<>>(t)
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
SEXP set_susceptibility_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {
  
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

[[cpp11::register]]
SEXP print_tool_cpp(SEXP t) {
  
  WrapTool(tptr)(t);
  tptr->print();
  return t;
  
}

// Function to get agent's viruses using get_agents_viruses()
[[cpp11::register]]
cpp11::writable::list get_agents_tools_cpp(SEXP model) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  
  cpp11::writable::list tools;
  
  for (auto & agent : ptr->get_agents())
    tools.push_back(
      cpp11::external_pointer< Tools<> >(
          new Tools<>(agent.get_tools())
      )
    );
  
  return tools;
  
}

[[cpp11::register]]
SEXP print_agent_tools_cpp(SEXP tools) {
  external_pointer<Tools<>> vptr(tools);
  vptr->print();
  return tools;
}

[[cpp11::register]]
SEXP set_distribution_tool_cpp(
  SEXP tool,
  SEXP model,
  SEXP tfun
  ) {
  
  WrapTool(toolptr)(tool);
  external_pointer<Model<>> mptr(model);
  external_pointer<ToolToAgentFun<>> tfunptr(tfun);
  
  toolptr->set_distribution(*tfunptr);
  
  return tool;
  
}

[[cpp11::register]]
SEXP distribute_tool_randomly_cpp(
  double prevalence,
  bool as_proportion
) {

  external_pointer<ToolToAgentFun<>> res(
      new ToolToAgentFun<>(
          distribute_tool_randomly(
            prevalence,
            as_proportion
          )
      )
  );
  
  return res;
  
}

[[cpp11::register]]
SEXP distribute_tool_to_set_cpp(
  integers agents_ids
) {

  // Converting integers to std::vector<size_t>
  std::vector<size_t> ids;
  for (auto & id : as_cpp<std::vector<int>>(agents_ids))
  {
    if (id < 0)
      stop("Agent's ID must be a positive integer.");
    ids.push_back(static_cast<size_t>(id));
  }


  external_pointer<ToolToAgentFun<>> res(
      new ToolToAgentFun<>(
          distribute_tool_to_set(ids)
      )
  );
  
  return res;
  
}



#undef WrapTool
