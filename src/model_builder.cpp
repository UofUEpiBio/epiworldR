#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"
#include "model-builder.h"

using namespace epiworld;

[[cpp11::register]]
SEXP Model_cpp() {
  return cpp11::external_pointer<Model<>>(
    new ModelBuilder()
  );
}

[[cpp11::register]]
SEXP add_state_cpp(SEXP model, std::string state_label, SEXP update_fun) {

  cpp11::external_pointer<Model<>> ptr(model);

  // Checking if the update function is NULL SEXP
  if (update_fun == R_NilValue) {
    ptr->add_state(state_label, nullptr);
    return model;
  }

  cpp11::external_pointer<UpdateFun<>> funptr(update_fun);
  ptr->add_state(state_label, *(funptr.get()));

  return model;

}


[[cpp11::register]]
SEXP update_fun_susceptible_cpp(cpp11::integers exclude) {

  // The empty exclusion list has the same semantics as the default sampler.
  // Retain its native function pointer so transmission dispatch uses the same
  // optimized path as models constructed directly in C++.
  if (exclude.size() == 0)
    return cpp11::external_pointer<UpdateFun<>>(
      new UpdateFun<>(default_update_susceptible<>)
    );

  std::vector< epiworld_fast_uint > exclude_vec;
  for (auto i : exclude) {
    if (i < 0)
      throw std::logic_error(
        "Trying to exclude a negative state: " + std::to_string(i) + "."
      );

    exclude_vec.push_back(static_cast<epiworld_fast_uint>(i));
  }

  return cpp11::external_pointer<UpdateFun<>>(
    new UpdateFun<>(
      std::move(
        sampler::make_update_susceptible<>(exclude_vec)
      )
    )
  );

}

[[cpp11::register]]
SEXP update_fun_infected_cpp() {
  return cpp11::external_pointer<UpdateFun<>>(
    new UpdateFun<>(default_update_exposed<>)
  );
}

[[cpp11::register]]
SEXP update_fun_incubation_cpp(int target_state) {
  if (target_state < 0)
    throw std::invalid_argument("target_state must be non-negative.");
  return cpp11::external_pointer<UpdateFun<>>(
    new UpdateFun<>([target_state](Agent<> * p, Model<> * m) {
      if (static_cast<size_t>(target_state) >= m->get_states().size())
        throw std::out_of_range("target_state is outside the model's states.");
      if (!p->get_virus())
        throw std::logic_error("The incubation update requires a virus.");
      if (m->runif() < 1.0 / p->get_virus()->get_incubation(m))
        p->change_state(*m, target_state);
    })
  );
}

[[cpp11::register]]
SEXP update_fun_recovery_cpp(std::string param_name) {
  return cpp11::external_pointer<UpdateFun<>>(
    new UpdateFun<>([rate = ParamRef(std::move(param_name))](Agent<> * p, Model<> * m) {
      if (!p->get_virus())
        throw std::logic_error("The recovery update requires a virus.");
      if (m->runif() < rate(*m))
        p->rm_virus(*m);
    })
  );
}

// Compatibility adapter for the existing model-specific initial proportions.
// The factories are shared core utilities, not model class instantiations.
[[cpp11::register]]
SEXP set_initial_states_recipe_cpp(SEXP model, std::string kind) {
  auto * ptr = dynamic_cast<ModelBuilder *>(
    cpp11::external_pointer<Model<>>(model).get()
  );
  if (!ptr)
    throw std::invalid_argument("Expected a model created with Model().");
  if (kind == "SIR")
    ptr->initializer_factory = epimodels::create_init_function_sir<int>;
  else if (kind == "SIRD")
    ptr->initializer_factory = epimodels::create_init_function_sird<int>;
  else if (kind == "SEIR")
    ptr->initializer_factory = epimodels::create_init_function_seir<int>;
  else if (kind == "SEIRD")
    ptr->initializer_factory = epimodels::create_init_function_seird<int>;
  else
    throw std::invalid_argument("Unknown initial-state recipe.");
  return model;
}

[[cpp11::register]]
SEXP update_fun_rate_cpp(
    cpp11::strings param_names,
    cpp11::integers target_states
) {

  std::vector< std::string > param_names_vec(param_names.begin(), param_names.end());
  std::vector< epiworld_fast_uint > target_states_vec;
  for (auto i : target_states) {
    if (i < 0)
      throw std::logic_error(
        "Trying to target a negative state: " + std::to_string(i) + "."
      );

    target_states_vec.push_back(static_cast<epiworld_fast_uint>(i));
  }

  return cpp11::external_pointer<UpdateFun<>>(
    new UpdateFun<>(
      std::move(new_state_update_transition<>(
        param_names_vec, target_states_vec
      ))
    )
  );

}

[[cpp11::register]]
SEXP set_state_function_cpp(SEXP model, std::string state_label, SEXP update_fun) {

  cpp11::external_pointer<Model<>> ptr(model);

  // Checking if the update function is NULL SEXP
  if (update_fun == R_NilValue) {
    ptr->set_state_function(state_label, nullptr);
    return model;
  }

  cpp11::external_pointer<UpdateFun<>> funptr(update_fun);
  ptr->set_state_function(state_label, *(funptr.get()));

  return model;

}
