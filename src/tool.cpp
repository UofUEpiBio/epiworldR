
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"
#include "assign.h"

using namespace cpp11;
using namespace epiworld;


// Model definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/models

#define WrapTool(a) cpp11::external_pointer<epiworld::Tool> (a)

[[cpp11::register]]
SEXP tool_cpp(
    std::string name,
    double susceptibility_reduction,
    double transmission_reduction,
    double recovery_enhancer,
    double death_reduction
    ) {

  WrapTool(tool)(new epiworld::Tool(name));

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
int add_tool_cpp(SEXP m, SEXP t, SEXP assign) {

  // A tool's placement takes no labels
  cpp11::external_pointer<epiworld::Model>(m)->add_tool(
    *cpp11::external_pointer<epiworld::Tool>(t),
    assign_from_r(assign, nullptr)
  );

  return 0;
}

[[cpp11::register]]
SEXP assign_tool_cpp(SEXP m, int tool_pos, SEXP assign) {

  cpp11::external_pointer<epiworld::Model>(m)->assign_tool(
    static_cast< size_t >(tool_pos), assign_from_r(assign, nullptr)
  );

  return m;

}

[[cpp11::register]]
SEXP rm_tool_cpp(SEXP m, size_t tool_pos) {
  cpp11::external_pointer<epiworld::Model>(m)->rm_tool(tool_pos);
  return m;
}

// Tool function ---------------------------------------------------------------
[[cpp11::register]]
SEXP tool_fun_logit_cpp(
    integers vars,
    doubles coefs,
    SEXP model
) {

  external_pointer<Model> mptr(model);

  external_pointer<ToolFun> res(
      new ToolFun(
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
  external_pointer<Model> mptr(model);

  toolptr->set_transmission_reduction(param);

  return tool;

}

[[cpp11::register]]
SEXP set_transmission_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {

  WrapTool(toolptr)(tool);
  external_pointer<Model> mptr(model);
  external_pointer<ToolFun> tfunptr(tfun);

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
  external_pointer<Model> mptr(model);

  toolptr->set_susceptibility_reduction(param);

  return tool;

}

[[cpp11::register]]
SEXP set_susceptibility_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {

  WrapTool(toolptr)(tool);
  external_pointer<Model> mptr(model);
  external_pointer<ToolFun> tfunptr(tfun);

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
  external_pointer<Model> mptr(model);

  toolptr->set_recovery_enhancer(
      param
  );

  return tool;

}

[[cpp11::register]]
SEXP set_recovery_enhancer_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {

  WrapTool(toolptr)(tool);
  external_pointer<Model> mptr(model);
  external_pointer<ToolFun> tfunptr(tfun);

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
  external_pointer<Model> mptr(model);

  toolptr->set_death_reduction(
      param
  );

  return tool;

}

[[cpp11::register]]
SEXP set_death_reduction_fun_cpp(SEXP tool, SEXP model, SEXP tfun) {

  WrapTool(toolptr)(tool);
  external_pointer<Model> mptr(model);
  external_pointer<ToolFun> tfunptr(tfun);

  toolptr->set_death_reduction_fun(*tfunptr);

  return tool;

}


[[cpp11::register]]
std::string get_name_tool_cpp(SEXP tool) {
  return external_pointer<Tool>(tool)->get_name();
}

[[cpp11::register]]
SEXP set_name_tool_cpp(SEXP tool, std::string name) {
  external_pointer<Tool>(tool)->set_name(name);
  return tool;
}

// Virus targets ---------------------------------------------------------------
[[cpp11::register]]
SEXP add_target_tool_cpp(SEXP tool, int lineage_id) {

  WrapTool(toolptr)(tool);
  toolptr->add_target(lineage_id);
  return tool;

}

[[cpp11::register]]
SEXP set_targets_tool_cpp(SEXP tool, integers lineage_ids) {

  WrapTool(toolptr)(tool);
  toolptr->set_targets(as_cpp<std::vector<int>>(lineage_ids));
  return tool;

}

[[cpp11::register]]
SEXP clear_targets_tool_cpp(SEXP tool) {

  WrapTool(toolptr)(tool);
  toolptr->clear_targets();
  return tool;

}

[[cpp11::register]]
cpp11::writable::integers get_targets_tool_cpp(SEXP tool) {

  WrapTool(toolptr)(tool);

  cpp11::writable::integers res;
  for (const auto & id : toolptr->get_targets())
    res.push_back(id);

  return res;

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

  cpp11::external_pointer<Model> ptr(model);

  cpp11::writable::list tools;

  for (auto & agent : ptr->get_agents())
    tools.push_back(
      cpp11::external_pointer< std::vector<ToolPtr> >(
          new std::vector<ToolPtr>(agent.get_tools())
      )
    );

  return tools;

}

[[cpp11::register]]
SEXP print_agent_tools_cpp(SEXP tools) {
  external_pointer<std::vector<ToolPtr> > vptr(tools);

  for (auto & tool : *vptr)
     tool->print();

  return tools;
}

#undef WrapTool

[[cpp11::register]]
SEXP clone_tool_cpp(SEXP tool) {

  // clone_ptr() keeps the tool's type (e.g., ToolVaccine)
  auto cloned = external_pointer<Tool>(tool)->clone_ptr();

  // The copy is a new tool: the model gives it its own id when it is added
  cloned->set_id(-99);

  external_pointer<Tool> ptr(cloned.release());

  return ptr;

}
