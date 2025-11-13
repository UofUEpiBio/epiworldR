#' Susceptible Exposed Infected Removed model (SEIR connected)
#'
#' The SEIR connected model implements a model where all agents are connected.
#' This is equivalent to a compartmental model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SEIR_model)).
#'
#' @param name String. Name of the virus.
#' @param n Integer greater than zero. Population size.
#' @param prevalence Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of
#' transmission.
#' @param incubation_days Numeric scalar greater than 0. Average number of
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery_rate.
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @concept general-models
#' @aliases epiworld_seirconn
#' @section Model diagram:
#' ![](seirconnected.png "SEIR Connected Diagram")
#' @returns
#' - The `ModelSEIRCONN`function returns a model of class [epiworld_model].
#' @examples
#' # An example with COVID-19
#' model_seirconn <- ModelSEIRCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01,
#'   n                   = 10000,
#'   contact_rate        = 2,
#'   incubation_days     = 7,
#'   transmission_rate   = 0.5,
#'   recovery_rate       = 0.3
#' )
#'
#' # Running and printing
#' run(model_seirconn, ndays = 100, seed = 1912)
#' model_seirconn
#'
#' plot(model_seirconn)
#'
#' # Adding the flu
#' flu <- virus("Flu", .9, 1 / 7, prevalence = 0.001, as_proportion = TRUE)
#' add_virus(model_seirconn, flu)
#'
#' #' # Running and printing
#' run(model_seirconn, ndays = 100, seed = 1912)
#' model_seirconn
#'
#' plot(model_seirconn)
#' @seealso epiworld-methods
ModelSEIRCONN <- function(
  name,
  n,
  prevalence,
  contact_rate,
  transmission_rate,
  incubation_days,
  recovery_rate
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(contact_rate)
  stopifnot_double(transmission_rate)
  stopifnot_double(incubation_days)
  stopifnot_double(recovery_rate)

  structure(
    ModelSEIRCONN_cpp(name, n, prevalence, contact_rate,
      transmission_rate, incubation_days, recovery_rate),
    class = c("epiworld_seirconn", "epiworld_model")
  )

}
