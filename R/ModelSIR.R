#' SIR model
#' @param name String. Name of the virus.
#'
#' Susceptible-Infected-Recovered model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIR_model)).
#'
#' @param name String. Name of the virus

#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @export
#' @family Models
#' @details
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @aliases epiworld_sir
#' @returns
#' - The `ModelSIR` function returns a model of class [epiworld_model].
#' @examples
#' model_sir <- ModelSIR(name = "COVID-19", prevalence = 0.01,
#'   transmission_rate = 0.9, recovery_rate = 0.1)
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_sir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
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
  prevalence,
  transmission_rate,
  recovery_rate
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_rate)
  stopifnot_double(recovery_rate)

  structure(
    ModelSIR_cpp(name, prevalence, transmission_rate, recovery_rate),
    class = c("epiworld_sir", "epiworld_model")
  )

}
