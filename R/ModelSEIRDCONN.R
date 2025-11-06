#' Susceptible Exposed Infected Removed Deceased model (SEIRD connected)
#'
#' The SEIRD connected model implements a model where all agents are connected.
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
#' @param death_rate Numeric scalar between 0 and 1. Probability of death.
#' @param n Number of individuals in the population.
#' @export
#' @details
#' The [initial_states] function allows the user to set the initial state of the
#' model. The user must provide a vector of proportions indicating the following
#' values: (1) Proportion of exposed agents who are infected, (2)
#' proportion of non-infected agents already removed, and (3) proportion of
#' non-ifected agents already deceased.
#' @section Model diagram:
#'
#' ![](diagrams/seirdconn.png "SEIRD Connected Diagram")
#' @family Models
#' @aliases epiworld_seirdconn
#' @returns
#' - The `ModelSEIRDCONN`function returns a model of class [epiworld_model].
#' @examples
#' # An example with COVID-19
#' model_seirdconn <- ModelSEIRDCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01,
#'   n                   = 10000,
#'   contact_rate        = 2,
#'   incubation_days     = 7,
#'   transmission_rate   = 0.5,
#'   recovery_rate       = 0.3,
#'   death_rate          = 0.01
#' )
#'
#' # Running and printing
#' run(model_seirdconn, ndays = 100, seed = 1912)
#' model_seirdconn
#'
#' plot(model_seirdconn)
#'
#' # Adding the flu
#' flu <- virus(
#'   "Flu", prob_infecting = .3, recovery_rate = 1 / 7,
#'   prob_death = 0.001,
#'   prevalence = 0.001, as_proportion = TRUE
#' )
#' add_virus(model = model_seirdconn, virus = flu)
#'
#' #' # Running and printing
#' run(model_seirdconn, ndays = 100, seed = 1912)
#' model_seirdconn
#'
#' plot(model_seirdconn)
#' @seealso epiworld-methods
ModelSEIRDCONN <- function(
    name,
    n,
    prevalence,
    contact_rate,
    transmission_rate,
    incubation_days,
    recovery_rate,
    death_rate
    ) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(contact_rate)
  stopifnot_double(transmission_rate)
  stopifnot_double(incubation_days)
  stopifnot_double(recovery_rate)
  stopifnot_double(death_rate)

  structure(
    ModelSEIRDCONN_cpp(name, n, prevalence, contact_rate,
      transmission_rate, incubation_days, recovery_rate,
      death_rate),
    class = c("epiworld_seirdconn", "epiworld_model")
  )

}
