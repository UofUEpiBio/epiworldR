#' Model building functions
#'
#' Functions to build models from scratch (or to modify existing models).
#' @name model_builder
#' @export
#' @examples
#' # Create a new model
#' model <- new_model()
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
#' virus_set_state(flu, 1, 2, 2)
#'
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
#' - The function `new_model()` returns a new model object of class
#' [epiworld_model].
new_model <- function() {

  new_model_cpp() |>
    structure(class = "epiworld_model")
}

#' @export
#' @param model An object of class [epiworld_model].
#' @param state_name A string with the name of the state to be added.
#' @param update_fun An object of class `epiworld_update_fun` with
#' the update function to be used for the new state. A `NULL` value
#' can be used if the state does not update (e.g., a "dead" state).
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
