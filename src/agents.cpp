
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"
#include "cpp11/list.hpp"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP get_agents_cpp(
  SEXP model
) {
  
  // Making some room
  
  cpp11::external_pointer<Model<>> ptr(model);
  cpp11::external_pointer<std::vector< Agent<> >> agents(
      &ptr->get_agents(),
      false
      );
  
  return agents;
  
}

[[cpp11::register]]
SEXP get_agent_cpp(
  SEXP agents,
  size_t i
) {
  
  cpp11::external_pointer<std::vector<Agent<>>> ptr(agents);
  
  if (i >= ptr->size())
    stop("The agent index %lu is out of range.\n", i);
  
  return cpp11::external_pointer< Agent<> >(
      new Agent<>(ptr->operator[](i))
    );

}

[[cpp11::register]]
SEXP print_agent_cpp(
    SEXP agent,
    SEXP model,
    bool compressed
) {
  
  cpp11::external_pointer<Agent<>> ptr(agent);
  cpp11::external_pointer<Model<>> ptr_model(model);
  
  ptr->print(&(*ptr_model), compressed);
  
  return agent;
  
}

[[cpp11::register]]
int get_state_agent_cpp(SEXP agent) {
  return cpp11::external_pointer<Agent<>>(agent)->get_state();
}

// Function to get agent's states using get_agents_states()
[[cpp11::register]]
std::vector<std::string> get_agents_states_cpp(SEXP model) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  
  std::vector<std::string> states;
  states.reserve(ptr->size());
  
  auto states_uint = ptr->get_agents_states();
  
  // Getting the model's states
  auto model_states = ptr->get_states();

  // Copying the data to states
  for (auto i : states_uint)
    states.push_back(model_states[i]);
  
  return states;
  
}
