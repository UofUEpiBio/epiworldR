#' SIRD model
#' @param name String. Name of the virus.
#'
#' Susceptible-Infected-Recovered-Deceased model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1166587407#The_SIRD_model)).
#'
#' @param name String. Name of the virus

#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @param death_rate Numeric scalar between 0 and 1. Rate of death from virus.
#' @export
#' @details
#' The [initial_states] function allows the user to set the initial state of the
#' model. The user must provide a vector of proportions indicating the following
#' values: (1) proportion of non-infected agents already removed, and (2) proportion of
#' non-ifected agents already deceased.
#' @family Models
#' @aliases epiworld_sird
#' @returns
#' - The `ModelSIRD` function returns a model of class [epiworld_model].
#' @examples
#' model_sird <- ModelSIRD(
#'   name              = "COVID-19",
#'   prevalence        = 0.01,
#'   transmission_rate = 0.9,
#'   recovery_rate     = 0.1,
#'   death_rate        = 0.01
#' )
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_sird,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
#'
#' # Running and printing
#' run(model_sird, ndays = 100, seed = 1912)
#' model_sird
#'
#' # Plotting
#' plot(model_sird)
#' @seealso epiworld-methods
ModelSIRD <- function(
    name, prevalence, transmission_rate, recovery_rate, death_rate
    ) {
  # Check inputs
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_rate)
  stopifnot_double(recovery_rate)
  stopifnot_double(death_rate)

  structure(
    ModelSIRD_cpp(name, prevalence, transmission_rate, recovery_rate, death_rate),
    class = c("epiworld_sird", "epiworld_model")
  )

}
