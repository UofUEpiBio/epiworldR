
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/matrix.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;

[[cpp11::register]]
SEXP get_entities_cpp(
  SEXP model
) {
  
  // Making some room
  
  cpp11::external_pointer<Model<>> ptr(model);

  cpp11::external_pointer<std::vector< Entity<> >> entities(
      &ptr->get_entities(),
      false
      );
  
  return entities;
  
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
  std::string name
) {
  
  cpp11::external_pointer<Entity<>> ptr(
    new Entity<>(name)
  );
  
  return ptr;
  
}

[[cpp11::register]]
int entity_size_cpp(SEXP entity) {
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
int entity_rm_agent_cpp(
  SEXP entity,
  int idx
  ) {
  
  cpp11::external_pointer<Entity<>> ptr(entity);
  
  ptr->rm_agent(static_cast<size_t>(idx)); 
  
  return 0;
}

[[cpp11::register]]
std::string entity_get_name_cpp(SEXP entity) {
  return cpp11::external_pointer<Entity<>>(entity)->get_name();
}

[[cpp11::register]]
int add_entity_cpp(SEXP model, SEXP entity) {
  
  cpp11::external_pointer<Model<>> ptr_model(model);
  cpp11::external_pointer<Entity<>> ptr_entity(entity);
  
  ptr_model->add_entity(*ptr_entity);
  
  return 0;
}

[[cpp11::register]]
int add_entity_n_cpp(SEXP model, SEXP entity, int n) {
  
  cpp11::external_pointer<Model<>> ptr_model(model);
  cpp11::external_pointer<Entity<>> ptr_entity(entity);
  
  ptr_model->add_entity_n(*ptr_entity, static_cast<size_t>(n));
  
  return 0;
}

// [[cpp11::register]]
// int rm_entity_cpp(SEXP model, int entity_pos) {
  
//   cpp11::external_pointer<Model<>> ptr_model(model);
  
//   ptr_model->rm_entity(
//     static_cast<size_t>(entity_pos)
//   );
  
//   return 0;
// }

// [[cpp11::register]]
// int entity_set_name_cpp(SEXP entity, std::string name) {
//   cpp11::external_pointer<Entity<>>(entity)->set_name(name);
//   return 0;
// }

