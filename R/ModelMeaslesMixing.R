#' Measles model with mixing
#'
#' `ModelMeaslesMixing` creates a measles epidemiological model with mixing
#' between different population groups. The model includes vaccination,
#' quarantine, isolation, and contact tracing mechanisms.
#'
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_matrix A row-stochastic matrix of mixing proportions between
#' population groups.
#' @param vax_reduction_recovery_rate Double. Vaccine reduction in recovery
#' rate (default: 0.5).
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of
#' transmission (default: 0.9).
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param prop_vaccinated Double. Proportion of population that is vaccinated.
#' @param vax_efficacy Double. Vaccine efficacy rate (default: 0.99).
#' @param quarantine_period Integer. Number of days for quarantine
#' (default: 21).
#' @param quarantine_willingness Double. Proportion of agents willing to
#' quarantine (default: 1).
#' @param isolation_willingness Double. Proportion of agents willing to isolate
#' (default: 1).
#' @param isolation_period Integer. Number of days for isolation (default: 4).
#' @param incubation_period Double. Duration of incubation period (default: 12).
#' @param prodromal_period Double. Duration of prodromal period (default: 4).
#' @param rash_period Double. Duration of rash period (default: 3).
#' @param hospitalization_rate Double. Rate of hospitalization (default: 0.2).
#' @param hospitalization_period Double. Period of hospitalization (default: 7).
#' @param days_undetected Double. Number of days an infection goes undetected
#' (default: 2).
#' @param contact_tracing_success_rate Double. Probability of successful
#' contact tracing (default: 1.0).
#' @param contact_tracing_days_prior Integer. Number of days prior to the onset
#' of the infection for which contact tracing is effective (default: 4).
#' @export
#' @family Models
#' @family measles models
#' @concept measles-models
#' @details
#' The `contact_matrix` is a matrix of contact rates between entities. The
#' matrix should be of size `n x n`, where `n` is the number of entities.
#' This is a row-stochastic matrix, i.e., the sum of each row should be 1.
#'
#' The model includes three distinct phases of measles infection: incubation,
#' prodromal, and rash periods. Vaccination provides protection against
#' infection and may reduce recovery time.
#'
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#'
#' The default value for the contact rate is an approximation to the disease's
#' basic reproduction number (R0), but it is not 100% accurate. A more accurate
#' way to se the contact rate is available, and will be distributed in the
#' future.
#' @section Model diagram:
#' ![](measlesmixing.png "Measles Mixing Diagram")
#' @returns
#' - The `ModelMeaslesMixing` function returns a model of classes
#' [epiworld_model] and [epiworld_measlesmixing].
#' @aliases epiworld_measlesmixing
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
#' measles_model <- ModelMeaslesMixing(
#'   n                        = N,
#'   prevalence               = 1 / N,
#'   contact_rate             = 15,
#'   transmission_rate        = 0.9,
#'   vax_efficacy             = 0.97,
#'   vax_reduction_recovery_rate = 0.8,
#'   incubation_period        = 10,
#'   prodromal_period         = 3,
#'   rash_period              = 7,
#'   contact_matrix           = cmatrix,
#'   hospitalization_rate     = 0.1,
#'   hospitalization_period   = 10,
#'   days_undetected          = 2,
#'   quarantine_period        = 14,
#'   quarantine_willingness   = 0.9,
#'   isolation_willingness    = 0.8,
#'   isolation_period         = 10,
#'   prop_vaccinated          = 0.95,
#'   contact_tracing_success_rate = 0.8,
#'   contact_tracing_days_prior = 4
#' )
#'
#' # Adding the entities to the model
#' measles_model |>
#'   add_entity(e1) |>
#'   add_entity(e2) |>
#'   add_entity(e3)
#'
#' set.seed(331)
#' run(measles_model, ndays = 100)
#' summary(measles_model)
#'
#' @seealso epiworld-methods
ModelMeaslesMixing <- function(
  n,
  prevalence,
  contact_matrix,
  vax_reduction_recovery_rate = .5,
  transmission_rate = .9,
  contact_rate = 15 / transmission_rate / prodromal_period,
  prop_vaccinated,
  vax_efficacy = .99,
  quarantine_period = 21,
  quarantine_willingness = 1,
  isolation_willingness = 1,
  isolation_period = 4,
  incubation_period = 12,
  prodromal_period = 4,
  rash_period = 3,
  hospitalization_rate = 0.2,
  hospitalization_period = 7,
  days_undetected = 2,
  contact_tracing_success_rate = 1.0,
  contact_tracing_days_prior = 4
) {
  # Check input parameters
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(contact_rate)
  stopifnot_double(transmission_rate)
  stopifnot_double(vax_efficacy, lb = 0, ub = 1)
  stopifnot_double(vax_reduction_recovery_rate)
  stopifnot_double(incubation_period)
  stopifnot_double(prodromal_period)
  stopifnot_double(rash_period)
  stopifany_na(contact_matrix)
  stopifnot_double(hospitalization_rate)
  stopifnot_double(hospitalization_period)
  stopifnot_double(days_undetected)
  stopifnot_int(quarantine_period)
  stopifnot_double(quarantine_willingness, lb = 0, ub = 1)
  stopifnot_double(isolation_willingness, lb = 0, ub = 1)
  stopifnot_int(isolation_period)
  stopifnot_double(prop_vaccinated, lb = 0, ub = 1)
  stopifnot_double(contact_tracing_success_rate, lb = 0, ub = 1)
  stopifnot_int(contact_tracing_days_prior, lb = 0)

  structure(
    ModelMeaslesMixing_cpp(
      n = n,
      prevalence = prevalence,
      contact_rate = contact_rate,
      transmission_rate = transmission_rate,
      vax_efficacy = vax_efficacy,
      vax_reduction_recovery_rate = vax_reduction_recovery_rate,
      incubation_period = incubation_period,
      prodromal_period = prodromal_period,
      rash_period = rash_period,
      contact_matrix = as.vector(contact_matrix),
      hospitalization_rate = hospitalization_rate,
      hospitalization_period = hospitalization_period,
      days_undetected = days_undetected,
      quarantine_period = quarantine_period,
      quarantine_willingness = quarantine_willingness,
      isolation_willingness = isolation_willingness,
      isolation_period = isolation_period,
      prop_vaccinated = prop_vaccinated,
      contact_tracing_success_rate = contact_tracing_success_rate,
      contact_tracing_days_prior = contact_tracing_days_prior
    ),
    class = c("epiworld_measlesmixing", "epiworld_model")
  )

}
