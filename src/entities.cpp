
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;

[[cpp11::register]]
SEXP get_entities_cpp(
  SEXP model
) {

  // Making some room

  cpp11::external_pointer<Model<>> ptr(model);

  cpp11::writable::list res;

  for (auto & i : ptr->get_entities()) {
    cpp11::external_pointer<Entity<>> entity(&i, false);
    res.push_back(entity);
  }

  return res;

}

[[cpp11::register]]
SEXP get_entity_cpp(
  SEXP entities,
  int idx
) {

  cpp11::external_pointer<std::vector< Entity<>> > ptr(entities);

  cpp11::external_pointer<Entity<>> entity(
    &ptr->at(static_cast<size_t>(idx)), false
  );

  return entity;

}


[[cpp11::register]]
SEXP entity_cpp(
  std::string name,
  double preval,
  bool as_proportion,
  bool to_unassigned
) {

  cpp11::external_pointer<Entity<>> ptr(
    new Entity<>(
      name,
      distribute_entity_randomly<>(
        preval,
        as_proportion,
        to_unassigned
      )
    )
  );

  return ptr;

}

[[cpp11::register]]
int get_entity_size_cpp(SEXP entity) {
  auto res = cpp11::external_pointer<Entity<>>(entity)->size();
  return static_cast<int>(res);
}

[[cpp11::register]]
int entity_add_agent_cpp(SEXP entity, SEXP agent, SEXP model) {

  cpp11::external_pointer<Entity<>> ptr(entity);
  cpp11::external_pointer<Agent<>> ptr_agent(agent);
  cpp11::external_pointer<Model<>> ptr_model(model);

  ptr->add_agent(&(*ptr_agent), &(*ptr_model));

  return 0;
}

[[cpp11::register]]
std::string get_entity_name_cpp(SEXP entity) {
  return cpp11::external_pointer<Entity<>>(entity)->get_name();
}

[[cpp11::register]]
int add_entity_cpp(
  SEXP model,
  SEXP entity
  ) {

  cpp11::external_pointer<Model<>> ptr_model(model);
  cpp11::external_pointer<Entity<>> ptr_entity(entity);

  ptr_model->add_entity(*ptr_entity);

  return 0;
}

[[cpp11::register]]
int rm_entity_cpp(SEXP model, int entity_pos) {

  cpp11::external_pointer<Model<>> ptr_model(model);

  ptr_model->rm_entity(
    static_cast<size_t>(entity_pos)
  );

  return 0;
}

[[cpp11::register]]
int load_agents_entities_ties_cpp(
  SEXP model,
  SEXP agents_ids,
  SEXP entities_ids
) {

  cpp11::external_pointer<Model<>> ptr_model(model);

  if (LENGTH(agents_ids) != LENGTH(entities_ids)) {
    cpp11::stop("agents_ids and entities_ids must have the same length");
  }

  ptr_model->load_agents_entities_ties(
    INTEGER(agents_ids),
    INTEGER(entities_ids),
    LENGTH(agents_ids)
  );

  return 0;

}

[[cpp11::register]]
cpp11::data_frame entity_get_agents_cpp(SEXP entity) {

  cpp11::external_pointer<Entity<>> ptr(entity);

  cpp11::writable::integers agent;
  cpp11::writable::integers entity_id;

  int id = static_cast<int>(ptr->get_id());
  for (const Agent<> & agent_i: ptr->get_agents()) {
    agent.push_back(static_cast<int>(agent_i.get_id()));
    entity_id.push_back(id);
  }

  return cpp11::writable::data_frame({
    "agent"_nm = agent,
    "entity"_nm = entity_id
  });

}

[[cpp11::register]]
int print_entity_cpp(SEXP entity) {
  cpp11::external_pointer<Entity<>>(entity)->print();
  return 0;
}

[[cpp11::register]]
SEXP set_distribution_entity_cpp(
  SEXP entity,
  SEXP fun
) {

  external_pointer<Entity<>> entity_ptr(entity);
  external_pointer<EntityToAgentFun<>> fun_ptr(fun);

  entity_ptr->set_distribution(*fun_ptr);

  return entity;

}

[[cpp11::register]]
SEXP distribute_entity_randomly_cpp(
  double prevalence,
  bool as_proportion,
  bool to_unassigned
) {

  external_pointer<EntityToAgentFun<>> res(
    new EntityToAgentFun<>(
      distribute_entity_randomly<>(
        prevalence,
        as_proportion,
        to_unassigned
      )
    )
  );

  return res;

}

[[cpp11::register]]
SEXP distribute_entity_to_set_cpp(
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

  external_pointer<EntityToAgentFun<>> res(
    new EntityToAgentFun<>(
      distribute_entity_to_set(ids)
    )
  );

  return res;

}

// [[cpp11::register]]
// int entity_set_name_cpp(SEXP entity, std::string name) {
//   cpp11::external_pointer<Entity<>>(entity)->set_name(name);
//   return 0;
// }
