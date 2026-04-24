#' Model building functions
#'
#' Functions to build models from scratch (or to modify existing models).
#'
#' @details
#' The model building functions allow users to create new models by adding
#' states, parameters, viruses, tools and other components. These functions are
#' useful for users who want to create custom models that are not included in
#' the package or to modify existing models.
#'
#' @name model_builder
#' @export
#' @examples
#' # Create a new model
#' model <- Model()
#'
#' # Adding recovery rate
#' add_param(model, "Rec Rate", 0.1)
#' add_param(model, "Trans Rate", 0.5)
#'
#' # Creating the update function for the susceptible state
#' update_fun_s <- update_fun_susceptible()
#'
#' # Adding the susceptible state to the model
#' add_state(model, "S", update_fun_s)
#'
#' # Creating the update function for the infected
#' update_fun_i <- update_fun_rate(
#'   param_names = "Rec Rate",
#'   target_states = 2L
#' )
#'
#' # Adding the infected state to the model
#' add_state(model, "I", update_fun_i)
#'
#' # Adding the recovered state to the model
#' add_state(model, "R", NULL)
#'
#' # Creating a virus
#' flu <- virus(
#'   "Flu", 0, .5, .2, .01,
#'   prevalence = 0, as_proportion = TRUE
#' )
#'
#' # We need to specify the effect that the virus
#' # has on agents when assigned, recovered, or removed.
#' # States are indexed starting at 0. In this model:
#' # 0 = S, 1 = I, 2 = R
#' virus_set_state(flu, 1, 2, 2)
#'
#' # We can set the transmission rate to be a function
#' # of the parameter "Trans Rate" using the
#' # set_prob_infecting_ptr function.
#' set_prob_infecting_ptr(flu, model, "Trans Rate")
#'
#' set_distribution_virus(
#'   flu,
#'   distribute_virus_randomly(1L, FALSE)
#' )
#'
#' # Adding the virus to the model
#' add_virus(model, flu)
#'
#' # Creating a SBM network
#' agents_smallworld(model, n = 1000, k = 20, d = FALSE, p = .01)
#'
#' # Running the model
#' run(model, ndays = 100, seed = 1912)
#'
#' summary(model)
#'
#' @return
#' - The function `Model()` returns a new model object of class
#' [epiworld_model].
#' @concept model-building-functions
Model <- function() {

  Model_cpp() |>
    structure(class = "epiworld_model")
}

#' @export
#' @param model An object of class [epiworld_model].
#' @param state_name A string with the name of the state to be added.
#' @param update_fun An object of class `epiworld_update_fun` with
#' the update function to be used for the new state. A `NULL` value
#' can be used if the state does not update (e.g., a "dead" state).
#' @return
#' - The function `add_state()` returns the modified model object
#' with the new state added. The function is called for its side
#' effects and returns the modified model invisibly.
#' @rdname model_builder
add_state <- function(
  model,
  state_name,
  update_fun = NULL
) {

  stopifnot_model(model)
  stopifnot_string(state_name)

  if (!is.null(update_fun))
    stopifnot_update_fun(update_fun)

  invisible(add_state_cpp(model, state_name, update_fun))

}

#' @export
#' @rdname model_builder
#' @param exclude An integer vector with the state indices to be
#' excluded from the infection process (see details).
#' @details
#' When using `update_fun_susceptible()`, the `exclude` argument can
#' be used to specify which agents carrying the virus should be excluded
#' from infecting susceptible agents. This can be useful in cases where
#' a virus may have a delayed effect on agents (e.g., an incubation period) or
#' when agents can recover but still carry the virus for some time.
#' @return
#' - The function `update_fun_susceptible()` returns an object of class
#' `epiworld_update_fun` that can be used as an update function for a
#' susceptible state.
#' @aliases epiworld_update_fun
update_fun_susceptible <- function(
  exclude = integer(0L)
) {

  stopifnot_int(exclude, lb = 0L)
  update_fun_susceptible_cpp(
    as.integer(exclude)
  ) |>
    structure(class = "epiworld_update_fun")

}

#' @export
#' @rdname model_builder
#' @param param_names A string vector with the name(s) of the parameter(s) to be
#' used in the update function.
#' @param target_states An integer vector with the state index(es) to which the
#' agent will transition when the update function is executed.
#' @details
#' The `update_fun_rate()` function creates an update function for a state that
#' updates at a constant rate (e.g., recovery, death). The `param_names`
#' argument specifies the name(s) of the parameter(s) that will be used in the
#' update function (e.g., "Recovery rate"). The `target_states` argument
#' specifies the state index(es) to which the agent will transition when the
#' update function is executed (e.g., 3 for a transition to a "Removed" state).
#' The function returns an object of class `epiworld_update_fun` that can be
#' used as an update function for a state in the model.
#'
#' Rates in the `update_fun_rate()` are daily probabilities of transitioning
#' to the target state(s). For example, if the recovery rate is 0.1, then there
#' is a 10% chance that an infected agent will transition to the "Removed"
#' state each day.
#' @return
#' - The function `update_fun_rate()` returns an object of class
#' `epiworld_update_fun` that can be used as an update function
#' for a state that updates at a constant rate (e.g., recovery, death).
update_fun_rate <- function(
  param_names,
  target_states
) {

  stopifnot_string(param_names)
  stopifnot_int(target_states, lb = 0L)
  update_fun_rate_cpp(
    param_names,
    as.integer(target_states)
  ) |>
    structure(class = "epiworld_update_fun")
}
