#' Building a model
#'
#' A model built from scratch: a population (`Model(n)`), its parameters
#' ([add_param()]), its contact layers ([add_layer()]) and its viruses
#' ([add_virus()]). The first virus gives the model its stages (see [virus()]).
#'
#' @param n Number of individuals in the population.
#' @param model An object of class `epiworld_model`.
#' @param assign A placement from the [assign] functions whose labels are
#' stages without the virus (e.g., `assign_shares(c(Recovered = 0.2))`).
#' @param transition Character scalar. A declared transition, by name (see
#' [add_external()]).
#' @param agent_id Integer scalar. The agent (starting at 0).
#'
#' @details
#' `set_initial()` places agents, at the start of each run, in stages without
#' the virus (e.g., recovered or vaccinated at the start). It draws from the
#' agents that are still in the first stage after the viruses are placed.
#'
#' A model without statuses has one, `"Free"`; the models that have statuses
#' (e.g., [ModelSEIRQuarantine()]) declare them in C++. `get_statuses()` lists
#' them.
#'
#' `fire()` applies a declared transition to an agent now (e.g., an external
#' transition from a global event). It returns whether the agent moved: it
#' does not if it is not in one of the transition's stages.
#'
#' @returns
#' - `Model()` returns a model, of class `epiworld_model`.
#' - `set_initial()` returns the model (invisibly).
#' - `get_statuses()` returns the labels of the model's statuses.
#' - `fire()` returns `TRUE` if the agent moved.
#' @examples
#' flu <- virus("Flu")
#' add_stage(flu, "Susceptible")
#' add_stage(flu, "Infected", infectious = TRUE)
#' add_stage(flu, "Recovered")
#' add_infection(flu, "Susceptible", "Infected", rate = "Transmission rate")
#' add_rate(flu, "Infected", "Recovered", rate = "Recovery rate")
#'
#' model <- Model(1000)
#' add_param(model, "Transmission rate", 0.3)
#' add_param(model, "Recovery rate", 0.1)
#' add_layer(model, layer_network("contacts", graph_smallworld(8, 0.1)))
#' add_virus(model, flu, assign_share(0.01))
#'
#' # A fifth of the agents start recovered
#' set_initial(model, assign_shares(c(Recovered = 0.2)))
#'
#' run(model, ndays = 50, seed = 1)
#' get_stages(model)
#' get_statuses(model)
#' @name model_builder
#' @export
#' @concept model-builder
Model <- function(n) {

  stopifnot_int(n)
  Model_cpp(n) |>
    structure(class = "epiworld_model")
}

#' @export
#' @rdname model_builder
set_initial <- function(model, assign) {

  stopifnot_model(model)
  stopifnot_assign(assign, null_ok = FALSE)
  set_initial_cpp(model, assign)
  if (uses_r(assign))
    attr(model, "single_threaded") <- TRUE
  invisible(model)

}

#' @export
#' @rdname model_builder
get_statuses <- function(model) {

  stopifnot_model(model)
  get_statuses_model_cpp(model)

}

#' @export
#' @rdname model_builder
fire <- function(model, transition, agent_id) {

  stopifnot_model(model)
  stopifnot_string(transition)
  stopifnot_int(agent_id, lb = 0L)
  fire_cpp(model, transition, as.integer(agent_id))

}
