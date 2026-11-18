#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"
#include "assign.h"
#include "natural-history.h"

using namespace epiworld;

[[cpp11::register]]
SEXP Model_cpp(unsigned int n) {
  return cpp11::external_pointer<Model>(
    new Model(n)
  );
}

// A disease's natural history (with several diseases, `disease` is required:
// get_natural_history() says so)
static const NaturalHistory & model_nh(const Model & model, int disease) {

  if (model.get_n_diseases() == 0u)
    cpp11::stop(
      "The model has no stages yet: add a virus with a natural history "
      "(add_virus())."
    );

  return *model.get_natural_history(disease);

}

[[cpp11::register]]
SEXP set_initial_cpp(SEXP model, SEXP assign, int disease) {

  cpp11::external_pointer<Model> ptr(model);
  const NaturalHistory & nh = model_nh(*ptr, disease);

  // The labels are the disease's stages
  ptr->set_initial(assign_from_r(assign, stage_resolver(nh)), disease);

  return model;

}

[[cpp11::register]]
cpp11::strings get_stages_model_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);
  model_nh(*ptr, disease);
  return cpp11::as_sexp(ptr->get_stages(disease));

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
bool fire_cpp(SEXP model, std::string transition, int agent_id, int disease) {

  cpp11::external_pointer<Model> ptr(model);
  return ptr->fire(
    ptr->transition_id(transition, disease),
    ptr->get_agent(static_cast< size_t >(agent_id))
  );

}

[[cpp11::register]]
std::string draw_natural_history_cpp(SEXP model) {
  return cpp11::external_pointer<Model>(model)->draw_natural_history();
}
