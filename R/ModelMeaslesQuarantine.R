#' Measles model with quarantine
#'
#' The `ModelMeaslesQuarantine` function implements a connected version
#' of a Measles model with quarantine.
#'
#' @param n Number of agents in the model.
#' @param prevalence Initial number of agents with the virus.
#' @param contact_rate Average number of contacts per step.
#' @param transmission_rate Probability of transmission.
#' @param vax_efficacy Probability of vaccine efficacy.
#' @param vax_improved_recovery Increase in recovery rate due to vaccination.
#' @param incubation_period Average number of incubation days.
#' @param prodromal_period Average number of prodromal days.
#' @param rash_period Average number of rash days.
#' @param days_undetected Average number of days undetected.
#' @param hospitalization_rate Probability of hospitalization.
#' @param hospitalization_period Average number of days in hospital.
#' @param prop_vaccinated Proportion of the population vaccinated.
#' @param quarantine_period Total duration of quarantine.
#' @param quarantine_willingness Probability of accepting quarantine (
#' see details).
#' @param isolation_period Average number of days in isolation.
#' @details
#' This model can be described as a SEIHR model with isolation and quarantine.
#' The infectious state is divided into prodromal and rash phases. Furthermore,
#' the quarantine state includes exposed, susceptible, prodromal, and recovered
#' states.
#'
#' The quarantine process is triggered any time that an agent with rash is
#' detected. The agent is then isolated and all agents who are unvaccinated are
#' quarantined. Isolated agents then may be moved out of the isolation in
#' `isolation_period` days.
#'
#' The basic reproductive number in Measles is estimated to be about 15.
#' By default, the contact rate of the model is set so that the R0 matches
#' 15.
#'
#' When `quarantine_period` is set to -1, the model assumes
#' there is no quarantine process. The same happens with `isolation_period`.
#' Since the quarantine process is triggered by an isolation, then
#' `isolation_period = -1` automatically sets `quarantine_period = -1`.
#' @export
#' @family Models
#' @aliases epiworld_measlesquarantine
#' @returns
#' - The `ModelMeaslesQuarantine` function returns a model of class [epiworld_model].
#' @examples
#' # An in a school with low vaccination
#' model_measles <- ModelMeaslesQuarantine(
#'   n = 500,
#'   prevalence = 1,
#'   prop_vaccinated = 0.70
#' )
#'
#' # Running and printing
#' run(model_measles, ndays = 100, seed = 1912)
#' model_measles
#'
#' plot(model_measles)
#'
#' @seealso epiworld-methods
#' @author
#' This model was built as a response to the US Measles outbreak in 2025.
#' This is a collaboration between the University of Utah (ForeSITE center
#' grant) and the Utah Department of Health and Human Services.
ModelMeaslesQuarantine <- function(
    n,
    prevalence = 1,
    contact_rate = 15 / transmission_rate / prodromal_period,
    transmission_rate = .9,
    vax_efficacy = .99,
    vax_improved_recovery = .5,
    incubation_period = 12,
    prodromal_period = 4,
    rash_period = 3,
    days_undetected = 2,
    hospitalization_rate = .2,
    hospitalization_period = 7,
    prop_vaccinated = 1 - 1 / 15,
    quarantine_period = 21,
    quarantine_willingness = 1,
    isolation_period = 4
    ) {
  # Check input parameters
  stopifnot_int(n, lb = 1)
  stopifnot_int(prevalence, lb = 0, ub = n)
  stopifnot_double(contact_rate, lb = 0)
  stopifnot_double(transmission_rate, lb = 0, ub = 1)
  stopifnot_double(vax_efficacy, lb = 0, ub = 1)
  stopifnot_double(vax_improved_recovery, lb = 0, ub = 1)
  stopifnot_double(incubation_period, lb = 0)
  stopifnot_double(prodromal_period, lb = 0)
  stopifnot_double(rash_period, lb = 0)
  stopifnot_double(days_undetected, lb = 0)
  stopifnot_double(hospitalization_rate, lb = 0, ub = 1)
  stopifnot_double(hospitalization_period, lb = 0)
  stopifnot_double(prop_vaccinated, lb = 0, ub = 1)
  stopifnot_int(quarantine_period, lb = -1)
  stopifnot_double(quarantine_willingness, lb = 0, ub = 1)
  stopifnot_int(isolation_period, lb = -1)

  structure(
    ModelMeaslesQuarantine_cpp(
      n,
      prevalence,
      contact_rate,
      transmission_rate,
      vax_efficacy,
      vax_improved_recovery,
      incubation_period,
      prodromal_period,
      rash_period,
      days_undetected,
      hospitalization_rate,
      hospitalization_period,
      prop_vaccinated,
      quarantine_period,
      quarantine_willingness
    ),
    class = c("epiworld_measlesquarantine", "epiworld_model")
  )

}
