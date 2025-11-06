#' SISD model
#'
#' Susceptible-Infected-Susceptible-Deceased model (SISD) ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIS_model))
#'
#' @param name String. Name of the virus.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery from virus.
#' @param death_rate Numeric scalar between 0 and 1. Rate of death from virus.
#' @export
#' @family Models
#' @returns
#' - The `ModelSISD` function returns a model of class [epiworld_model].
#' @aliases epiworld_sisd
#' @examples
#' model_sisd <- ModelSISD(
#'   name = "COVID-19",
#'   prevalence = 0.01,
#'   transmission_rate = 0.9,
#'   recovery_rate = 0.1,
#'   death_rate = 0.01
#' )
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_sisd,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
#'
#' # Running and printing
#' run(model_sisd, ndays = 100, seed = 1912)
#' model_sisd
#'
#' # Plotting
#' plot(model_sisd, main = "SISD Model")
#'
#' @seealso epiworld-methods
ModelSISD <- function(
  name,
  prevalence,
  transmission_rate,
  recovery_rate,
  death_rate) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(transmission_rate)
  stopifnot_double(recovery_rate)
  stopifnot_double(death_rate)

  structure(
    ModelSISD_cpp(name, prevalence, transmission_rate, recovery_rate, death_rate),
    class = c("epiworld_sisd", "epiworld_model")
  )

}
