#' Measles model with mixing
#'
#' `ModelMeaslesMixing` creates a measles epidemiological model with mixing
#' between different population groups. The model includes vaccination,
#' quarantine, isolation, and contact tracing mechanisms.
#'
#' @param vname String. Name of the virus
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of
#' transmission.
#' @param vax_efficacy Double. Vaccine efficacy rate.
#' @param vax_reduction_recovery_rate Double. Vaccine reduction in recovery rate.
#' @param incubation_period Double. Duration of incubation period.
#' @param prodromal_period Double. Duration of prodromal period.
#' @param rash_period Double. Duration of rash period.
#' @param contact_matrix Matrix of contact rates between individuals.
#' @param hospitalization_rate Double. Rate of hospitalization.
#' @param hospitalization_period Double. Period of hospitalization.
#' @param days_undetected Double. Number of days an infection goes undetected.
#' @param quarantine_period Integer. Number of days for quarantine.
#' @param quarantine_willingness Double. Proportion of agents willing to quarantine.
#' @param isolation_willingness Double. Proportion of agents willing to isolate.
#' @param isolation_period Integer. Number of days for isolation.
#' @param prop_vaccinated Double. Proportion of population that is vaccinated.
#' @param contact_tracing_success_rate Double. Probability of successful contact tracing.
#' @param contact_tracing_days_prior Integer. Number of days prior to the onset
#' of the infection for which contact tracing is effective.
#' @export
#' @family Models
#' @details
#' The `contact_matrix` is a matrix of contact rates between entities. The
#' matrix should be of size `n x n`, where `n` is the number of entities.
#' This is a row-stochastic matrix, i.e., the sum of each row should be 1.
#'
#' The model includes three distinct phases of measles infection: incubation,
#' prodromal, and rash periods. Vaccination provides protection against infection
#' and may reduce recovery time.
#'
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @returns
#' - The `ModelMeaslesMixing` function returns a model of class [epiworld_model].
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
#'   vname                    = "Measles",
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
    vname,
    n,
    prevalence,
    contact_rate,
    transmission_rate,
    vax_efficacy,
    vax_reduction_recovery_rate,
    incubation_period,
    prodromal_period,
    rash_period,
    contact_matrix,
    hospitalization_rate,
    hospitalization_period,
    days_undetected,
    quarantine_period,
    quarantine_willingness,
    isolation_willingness,
    isolation_period,
    prop_vaccinated,
    contact_tracing_success_rate = 1.0,
    contact_tracing_days_prior = 4
    ) {
  # Check input parameters
  stopifnot_string(vname)
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
      vname, n, prevalence, contact_rate,
      transmission_rate, vax_efficacy,
      vax_reduction_recovery_rate, incubation_period,
      prodromal_period, rash_period,
      as.vector(contact_matrix),
      hospitalization_rate, hospitalization_period,
      days_undetected, quarantine_period,
      quarantine_willingness,
      isolation_willingness,
      isolation_period,
      prop_vaccinated,
      contact_tracing_success_rate,
      contact_tracing_days_prior
    ),
    class = c("epiworld_measlesmixing", "epiworld_model")
  )

}
