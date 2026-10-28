#' Susceptible Exposed Infected Recovered model (SEIR)
#'
#' @param name String. Name of the virus.
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param latent_period Numeric scalar greater than 0. Average number of
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @export
#' @family Models
#' @concept general-models
#' @aliases epiworld_seir
#' @details
#' The [initial_states] function allows the user to set the initial state of the
#' model. The user must provide a vector of proportions indicating the following
#' values: (1) Proportion of non-infected agents who are removed, and (2)
#' Proportion of exposed agents to be set as infected.
#' @section Model diagram:
#' ![](seir.png "SEIR Diagram")
#' @returns
#' - The `ModelSEIR`function returns a model of class [epiworld_model].
#' @examples
#' model_seir <- ModelSEIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1, latent_period = 4)
#'
#' # Adding a small world population
#' add_layer(model_seir, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_seir, ndays = 100, seed = 1912)
#' model_seir
#'
#' plot(model_seir, main = "SEIR Model")
#' @seealso epiworld-methods
ModelSEIR <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  latent_period,
  recovery_rate
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_prob)
  stopifnot_double(latent_period)
  stopifnot_double(recovery_rate)

  structure(
    ModelSEIR_cpp(name, n, prevalence, transmission_prob, latent_period, recovery_rate),
    class = c("epiworld_seir", "epiworld_model")
  )

}
