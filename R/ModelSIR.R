#' SIR model
#' @param name String. Name of the virus.
#'
#' Susceptible-Infected-Recovered model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIR_model)).
#'
#' @param name String. Name of the virus

#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @export
#' @family Models
#' @concept general-models
#' @details
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @section Model diagram:
#' ![](sir.png "SIR Diagram")
#' @aliases epiworld_sir
#' @returns
#' - The `ModelSIR` function returns a model of class [epiworld_model].
#' @examples
#' model_sir <- ModelSIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1)
#'
#' # Adding a small world population
#' add_layer(model_sir, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_sir, ndays = 100, seed = 1912)
#' model_sir
#'
#' # Plotting
#' plot(model_sir)
#' @seealso epiworld-methods
ModelSIR <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  recovery_rate
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_prob)
  stopifnot_double(recovery_rate)

  structure(
    ModelSIR_cpp(name, n, prevalence, transmission_prob, recovery_rate),
    class = c("epiworld_sir", "epiworld_model")
  )

}
