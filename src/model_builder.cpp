#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"

using namespace epiworld;

[[cpp11::register]]
SEXP Model_cpp() {
  return cpp11::external_pointer<Model<>>(
    new Model<>()
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
