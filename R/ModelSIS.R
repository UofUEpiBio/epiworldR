#' SIS model
#'
#' Susceptible-Infected-Susceptible model (SIS) ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIS_model))
#'
#' @param name String. Name of the virus.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery from virus.
#' @export
#' @family Models
#' @concept general-models
#' @returns
#' - The `ModelSIS` function returns a model of class [epiworld_model].
#' @aliases epiworld_sis
#' @section Model diagram:
#' ![](diagrams/sis.png "SIS Diagram")
#' @examples
#' model_sis <- ModelSIS(name = "COVID-19", prevalence = 0.01,
#'   transmission_rate = 0.9, recovery_rate = 0.1)
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_sis,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
#'
#' # Running and printing
#' run(model_sis, ndays = 100, seed = 1912)
#' model_sis
#'
#' # Plotting
#' plot(model_sis, main = "SIS Model")
#'
#' @seealso epiworld-methods
ModelSIS <- function(
  name,
  prevalence,
  transmission_rate,
  recovery_rate) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_rate)
  stopifnot_double(recovery_rate)

  structure(
    ModelSIS_cpp(name, prevalence, transmission_rate, recovery_rate),
    class = c("epiworld_sis", "epiworld_model")
  )

}
