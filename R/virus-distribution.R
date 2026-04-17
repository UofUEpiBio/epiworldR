#' Virus distribution functions
#'
#' Distribution functions control how a virus is initially assigned across
#' agents in a model.
#'
#' @param virus An object of class `epiworld_virus`.
#' @param distfun An object of class `epiworld_virus_distfun`.
#' @param prevalence Numeric scalar. Prevalence of the virus. In the case of
#' `distribute_virus_to_entities()`, it is a vector of prevalences, one per
#' entity.
#' @param as_proportion Logical scalar. If `TRUE`, `prevalence` is interpreted
#' as a proportion of the total number of agents in the model.
#' @param agents_ids Integer vector. Indices of the agents that will receive the
#' virus.
#' @returns
#' - `set_distribution_virus()` does not return a value. It assigns a
#' distribution function to the specified virus.
#' - `distribute_virus_randomly()` returns a distribution function of class
#' `epiworld_virus_distfun`. When `agents_ids` is not empty, it distributes the
#' virus randomly within that set. Otherwise it uses all agents in the model.
#' - `distribute_virus_to_set()` returns a distribution function of class
#' `epiworld_virus_distfun`.
#' - `distribute_virus_set()` is a deprecated alias for
#' `distribute_virus_to_set()`.
#' - `distribute_virus_to_entities()` returns a distribution function of class
#' `epiworld_virus_distfun`.
#' @details
#' `set_distribution_virus()` assigns a distribution function to the specified
#' virus of class [epiworld_virus]. Distribution functions can be created with
#' [distribute_virus_randomly()], [distribute_virus_to_set()], and
#' [distribute_virus_to_entities()].
#'
#' `distribute_virus_randomly()` creates a distribution function that randomly
#' assigns the virus to a proportion of the population.
#'
#' `distribute_virus_to_set()` creates a distribution function that assigns the
#' virus to a fixed set of agents.
#'
#' `distribute_virus_to_entities()` creates a distribution function that
#' assigns the virus to a number of agents based on entity-level prevalence.
#' This is only useful for the mixing models.
#'
#' @examples
#' model <- ModelSEIRCONN(
#'   name              = "COVID-19",
#'   prevalence        = 0.01,
#'   n                 = 10000,
#'   contact_rate      = 4,
#'   incubation_days   = 7,
#'   transmission_rate = 0.5,
#'   recovery_rate     = 0.99
#' )
#'
#' delta <- virus(
#'   "Delta Variant", 0.3, TRUE, .5, .2, .01
#' )
#'
#' set_distribution_virus(
#'   delta,
#'   distribute_virus_randomly(100, as_proportion = FALSE)
#' )
#'
#' add_virus(model, delta)
#'
#' set_distribution_virus(
#'   delta,
#'   distribute_virus_to_set(1:10)
#' )
#'
#' @name virus-distribution
#' @concept virus-distribution
NULL

#' @export
#' @rdname virus-distribution
set_distribution_virus <- function(virus, distfun) {

  stopifnot_virus(virus)
  stopifnot_virus_distfun(distfun)
  invisible(set_distribution_virus_cpp(virus, distfun))

}

#' @export
#' @rdname virus-distribution
distribute_virus_randomly <- function(
  prevalence,
  as_proportion,
  agents_ids = integer(0)
) {

  structure(
    distribute_virus_randomly_cpp(
      as.double(prevalence),
      as.logical(as_proportion),
      as.integer(agents_ids)
    ),
    class = "epiworld_virus_distfun"
  )

}

#' @export
#' @rdname virus-distribution
distribute_virus_to_set <- function(agents_ids) {

  structure(
    distribute_virus_to_set_cpp(as.vector(agents_ids)),
    class = "epiworld_virus_distfun"
  )

}

#' @export
#' @rdname virus-distribution
distribute_virus_set <- function(agents_ids) {

  .Deprecated("distribute_virus_to_set")

}

#' @export
#' @rdname virus-distribution
distribute_virus_to_entities <- function(
  prevalence,
  as_proportion
) {

  structure(
    distribute_virus_to_entities_cpp(
      as.double(prevalence),
      as.logical(as_proportion)
    ),
    class = "epiworld_virus_distfun"
  )

}
