#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"
#include "assign.h"

using namespace epiworld;

[[cpp11::register]]
SEXP Model_cpp(unsigned int n) {
  return cpp11::external_pointer<Model>(
    new Model(n)
  );
}

// The model's natural history (its stages come from its first virus)
static const NaturalHistory & model_nh(const Model & model) {

  const NaturalHistory * nh = model.get_natural_history();
  if (nh == nullptr)
    cpp11::stop(
      "The model has no stages yet: add a virus with a natural history "
      "(add_virus())."
    );

  return *nh;

}

[[cpp11::register]]
SEXP set_initial_cpp(SEXP model, SEXP assign) {

  cpp11::external_pointer<Model> ptr(model);
  const NaturalHistory & nh = model_nh(*ptr);

  // The labels are the model's stages
  ptr->set_initial(
    assign_from_r(assign, [&nh](const std::string & label) {
      return nh.stage_id(label);
    })
  );

  return model;

}

[[cpp11::register]]
cpp11::strings get_stages_model_cpp(SEXP model) {

  cpp11::external_pointer<Model> ptr(model);
  return cpp11::as_sexp(model_nh(*ptr).get_stages());

}

[[cpp11::register]]
cpp11::strings get_statuses_model_cpp(SEXP model) {

  cpp11::external_pointer<Model> ptr(model);
  const Statuses * st = ptr->get_statuses();
  if (st == nullptr)
    return cpp11::as_sexp(std::vector< std::string >({"Free"}));

  return cpp11::as_sexp(st->get_stages());

}

[[cpp11::register]]
bool fire_cpp(SEXP model, std::string transition, int agent_id) {

  cpp11::external_pointer<Model> ptr(model);
  return ptr->fire(
    ptr->transition_id(transition),
    ptr->get_agent(static_cast< size_t >(agent_id))
  );

}

[[cpp11::register]]
std::string draw_natural_history_cpp(SEXP model) {
  return cpp11::external_pointer<Model>(model)->draw_natural_history();
}
