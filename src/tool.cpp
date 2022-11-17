
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace cpp11;


// Model definitions:
// https://github.com/UofUEpi/epiworld/tree/master/include/epiworld/models

#define WrapTool(a) \
  cpp11::external_pointer<epiworld::Tool<>> (a)

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
  
#undef WrapTool
