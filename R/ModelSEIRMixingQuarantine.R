#' Susceptible Exposed Infected Removed model (SEIR) with mixing and quarantine
#'
#' `ModelSEIRMixingQuarantine` creates a model of the SEIR type with mixing and
#' a quarantine mechanism. Agents who are infected can be quarantined or
#' isolated. Isolation happens after the agent has been detected as infected,
#' and agents who have been in contact with the detected person will me moved
#' to quarantined status.
#'
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of
#' transmission.
#' @param incubation_days Numeric scalar. Average number of days in the
#' incubation period.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery.
#' @param n Number of individuals in the population.
#' @param contact_matrix Matrix of contact rates between individuals.
#' @param hospitalization_rate Double. Rate of hospitalization.
#' @param hospitalization_period Double. Period of hospitalization.
#' @param days_undetected Double. Number of days an infection goes undetected.
#' @param quarantine_period Integer. Number of days for quarantine.
#' @param quarantine_willingness Double. Proportion of agents willing to quarantine.
#' @param isolation_willingness Double. Proportion of agents willing to isolate.
#' @param isolation_period Integer. Number of days for isolation.
#' @param contact_tracing_success_rate Double. Probability of successful contact tracing.
#' @param contact_tracing_days_prior Integer. Number of days prior to the onset
#' of the infection for which contact tracing is effective.
#' @export
#' @family Models
#' @concept general-models
#' @details
#' The `contact_matrix` is a matrix of contact rates between entities. The
#' matrix should be of size `n x n`, where `n` is the number of entities.
#' This is a row-stochastic matrix, i.e., the sum of each row should be 1.
#'
#' The quarantine and isolation processes can be turned off by specifying
#' `quarantine_period < 0` and `isolation_period < 0` respectively. In other
#' words, any negative value for these parameters will suppress the
#' corresponding process.
#'
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @section Model diagram:
#' ![](seirmixingquarantine.png "SEIR Mixing Quarantine Diagram")
#' @returns
#' - The `ModelSEIRMixingQuarantine` function returns a model of class [epiworld_model].
#' @aliases epiworld_seirmixingquarantine
#'
#' @examples
#'
#' # Start off creating three entities.
#' # Individuals will be distributed randomly between the three.
#' e1 <- entity("Population 1", 3e3, as_proportion = FALSE)
#' e2 <- entity("Population 2", 3e3, as_proportion = FALSE)
#' e3 <- entity("Population 3", 3e3, as_proportion = FALSE)
#'
#' # Row-stochastic matrix (rowsums 1)
#' cmatrix <- c(
#'   c(0.9, 0.05, 0.05),
#'   c(0.1, 0.8, 0.1),
#'   c(0.1, 0.2, 0.7)
#' ) |> matrix(byrow = TRUE, nrow = 3)
#'
#' N <- 9e3
#'
#' flu_model <- ModelSEIRMixingQuarantine(
#'   name                  = "Flu",
#'   n                     = N,
#'   prevalence            = 1 / N,
#'   contact_rate          = 20,
#'   transmission_rate     = 0.1,
#'   recovery_rate         = 1 / 7,
#'   incubation_days       = 7,
#'   contact_matrix        = cmatrix,
#'   hospitalization_rate  = 0.05,
#'   hospitalization_period = 7,
#'   days_undetected       = 3,
#'   quarantine_period     = 14,
#'   quarantine_willingness = 0.8,
#'   isolation_period      = 7,
#'   isolation_willingness = 0.5,
#'   contact_tracing_success_rate = 0.7,
#'   contact_tracing_days_prior = 3
#' )
#'
#' # Adding the entities to the model
#' flu_model |>
#'   add_entity(e1) |>
#'   add_entity(e2) |>
#'   add_entity(e3)
#'
#' set.seed(331)
#' run(flu_model, ndays = 100)
#' summary(flu_model)
#'
#' @seealso epiworld-methods
ModelSEIRMixingQuarantine <- function(
  name,
  n,
  prevalence,
  contact_rate,
  transmission_rate,
  incubation_days,
  recovery_rate,
  contact_matrix,
  hospitalization_rate,
  hospitalization_period,
  days_undetected,
  quarantine_period,
  quarantine_willingness,
  isolation_willingness,
  isolation_period,
  contact_tracing_success_rate,
  contact_tracing_days_prior
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n, lb = 1L)
  stopifnot_double(prevalence, lb = 0, ub = 1)
  stopifnot_double(contact_rate, lb = 0)
  stopifnot_double(transmission_rate, lb = 0, ub = 1)
  stopifnot_double(incubation_days, lb = 0)
  stopifnot_double(recovery_rate, lb = 0, ub = 1)
  stopifany_na(contact_matrix)
  stopifnot_double(hospitalization_rate, lb = 0, ub = 1)
  stopifnot_double(hospitalization_period, lb = 0)
  stopifnot_double(days_undetected)
  stopifnot_int(quarantine_period)
  stopifnot_double(quarantine_willingness, lb = 0, ub = 1)
  stopifnot_double(isolation_willingness, lb = 0, ub = 1)
  stopifnot_int(isolation_period)
  stopifnot_double(contact_tracing_success_rate, lb = 0, ub = 1)
  stopifnot_int(contact_tracing_days_prior, lb = 0)

  structure(
    ModelSEIRMixingQuarantine_cpp(
      name, n, prevalence, contact_rate,
      transmission_rate, incubation_days,
      recovery_rate, as.vector(contact_matrix),
      hospitalization_rate, hospitalization_period,
      days_undetected, quarantine_period,
      quarantine_willingness,
      isolation_willingness,
      isolation_period,
      contact_tracing_success_rate,
      contact_tracing_days_prior
    ),
    class = c("epiworld_seirmixingquarantine", "epiworld_model")
  )

}
