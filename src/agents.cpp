
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

  cpp11::external_pointer<Model> ptr(model);
  cpp11::external_pointer<std::vector< Agent >> agents(
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

  cpp11::external_pointer<std::vector<Agent>> ptr(agents);

  if (i >= ptr->size())
    stop("The agent index %lu is out of range.\n", i);

  // The model's own agent, not a copy: add_virus_agent() and add_tool_agent()
  // act on it (the R object keeps the model alive)
  return cpp11::external_pointer< Agent >(&ptr->operator[](i), false);

}

[[cpp11::register]]
SEXP print_agent_cpp(
    SEXP agent,
    SEXP model,
    bool compressed
) {

  cpp11::external_pointer<Agent>(agent)->
    print(
      *cpp11::external_pointer<Model>(model),
      compressed
    );
  return agent;

}

// The agent's stage and status, by label
[[cpp11::register]]
std::string get_stage_agent_cpp(SEXP model, SEXP agent, int disease) {

  cpp11::external_pointer<Model> ptr(model);
  const int s = ptr->get_stage(*cpp11::external_pointer<Agent>(agent), disease);
  return ptr->get_stages(disease)[static_cast<size_t>(s)];

}

[[cpp11::register]]
std::string get_status_agent_cpp(SEXP model, SEXP agent) {

  cpp11::external_pointer<Model> ptr(model);
  const Statuses * st = ptr->get_statuses();
  if (st == nullptr)
    return "Free";

  const int k = ptr->get_status(*cpp11::external_pointer<Agent>(agent));
  return st->get_stages()[static_cast<size_t>(k)];

}

// Every agent's stage (or status), by label
[[cpp11::register]]
std::vector<std::string> get_agents_stages_cpp(
  SEXP model,
  bool status,
  int disease
) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector<std::string> res;
  res.reserve(ptr->size());

  if (status)
  {
    const Statuses * st = ptr->get_statuses();
    for (auto & a : ptr->get_agents())
      res.push_back(st == nullptr ? "Free" :
        st->get_stages()[static_cast<size_t>(ptr->get_status(a))]);
    return res;
  }

  if (ptr->get_n_diseases() == 0u)
  {
    res.assign(ptr->size(), "");
    return res;
  }

  const auto & labels = ptr->get_stages(disease);
  for (int s : ptr->get_agents_stages(disease))
    res.push_back(labels[static_cast<size_t>(s)]);

  return res;

}

[[cpp11::register]]
SEXP add_virus_agent_cpp(SEXP agent, SEXP model, SEXP virus) {


  cpp11::external_pointer<Agent> ptr_agent(agent);
  cpp11::external_pointer<Model> ptr_model(model);
  cpp11::external_pointer<Virus> ptr_virus(virus);

  ptr_agent->set_virus(*ptr_model, *ptr_virus);

  return agent;

}

[[cpp11::register]]
SEXP add_tool_agent_cpp(SEXP agent, SEXP model, SEXP tool) {

  cpp11::external_pointer<Agent> ptr_agent(agent);
  cpp11::external_pointer<Model> ptr_model(model);
  cpp11::external_pointer<Tool> ptr_tool(tool);

  ptr_agent->add_tool(*ptr_model, *ptr_tool);

  return agent;

}


[[cpp11::register]]
bool has_virus_cpp(SEXP agent, SEXP virus) {

  cpp11::external_pointer<Agent> ptr_agent(agent);
  cpp11::external_pointer<Virus> ptr_virus(virus);

  return ptr_agent->has_virus(*ptr_virus);

}

[[cpp11::register]]
bool has_tool_cpp(SEXP agent, SEXP tool) {

  cpp11::external_pointer<Agent> ptr_agent(agent);
  cpp11::external_pointer<Tool> ptr_tool(tool);

  return ptr_agent->has_tool(*ptr_tool);

}
