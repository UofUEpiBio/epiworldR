#' Measles model with quarantine
#'
#' Implements a Susceptible-Exposed-Infectious-Hospitalized-Recovered (SEIHR)
#' model for Measles within a school. The model includes isolation of
#' detected cases and optional quarantine of unvaccinated individuals.
#'
#' @param n Number of agents in the model.
#' @param prevalence Initial number of agents with the virus.
#' @param contact_rate Average number of contacts per step. Default is set to
#' match the basic reproductive number (R0) of 15 (see details).
#' @param transmission_rate Probability of transmission.
#' @param vax_efficacy Probability of vaccine efficacy.
#' @param incubation_period Average number of incubation days.
#' @param prodromal_period Average number of prodromal days.
#' @param rash_period Average number of rash days.
#' @param days_undetected Average number of days undetected. Detected cases
#' are moved to isolation and trigger the quarantine process.
#' @param hospitalization_rate Probability of hospitalization.
#' @param hospitalization_period Average number of days in hospital.
#' @param prop_vaccinated Proportion of the population vaccinated.
#' @param quarantine_period Number of days an agent is in quarantine.
#' @param quarantine_willingness Probability of accepting quarantine (
#' see details).
#' @param isolation_period Number of days an agent is in isolation.
#' @param ... Further arguments (not used).
#' @details
#' This model can be described as a SEIHR model with isolation and quarantine.
#' The infectious state is divided into prodromal and rash phases. Furthermore,
#' the quarantine state includes exposed, susceptible, prodromal, and recovered
#' states.
#'
#' The model is a perfect mixing model, meaning that all agents are in contact
#' with each other. The model is designed to simulate the spread of Measles
#' within a school setting, where the population is assumed to be homogeneous.
#'
#' The quarantine process is triggered any time that an agent with rash is
#' detected. The agent is then isolated and all agents who are unvaccinated are
#' quarantined (if willing). Isolated agents then may be moved out of the
#' isolation in `isolation_period` days. The quarantine willingness parameter
#' sets the probability of accepting quarantine. If a quarantined agent develops
#' rash, they are moved to isolation, which triggers a new quarantine process.
#'
#' The basic reproductive number in Measles is estimated to be about 15.
#' By default, the contact rate of the model is set so that the R0 matches
#' 15.
#'
#' When `quarantine_period` is set to -1, the model assumes
#' there is no quarantine process. The same happens with `isolation_period`.
#' Since the quarantine process is triggered by an isolation, then
#' `isolation_period = -1` automatically sets `quarantine_period = -1`.
#'
#' @note
#' As of version 0.10.0, the parameter `vax_improved_recovery` has been removed
#' and is no longer used (it never had a side effect). Future versions may not
#' accept it.
#'
#' @references
#' Jones, Trahern W, and Katherine Baranowski. 2019. "Measles and Mumps: Old
#' Diseases, New Outbreaks."
#'
#' Liu, Fengchen, Wayne T A Enanoria, Jennifer Zipprich, Seth Blumberg,
#' Kathleen Harriman, Sarah F Ackley, William D Wheaton, Justine L
#' Allpress, and Travis C Porco. 2015. "The Role of Vaccination Coverage,
#' Individual Behaviors, and the Public Health Response in the Control of
#' Measles Epidemics: An Agent-Based Simulation for California." *BMC
#' Public Health* 15 (1): 447. \doi{10.1186/s12889-015-1766-6}.
#'
#' "Measles Disease Plan." 2019. Utah Department of Health and Human
#' Services. <https://epi.utah.gov/wp-content/uploads/Measles-disease-plan.pdf>.
#' @export
#' @family Models
#' @family measles models
#' @aliases epiworld_measlesquarantine
#' @returns
#' - The `ModelMeaslesQuarantine` function returns a model of classes [epiworld_model] and `epiworld_measlesquarantine`.
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
ModelMeaslesSchool <- function(
    n,
    prevalence = 1,
    contact_rate = 15 / transmission_rate / prodromal_period,
    transmission_rate = .9,
    vax_efficacy = .99,
    incubation_period = 12,
    prodromal_period = 4,
    rash_period = 3,
    days_undetected = 2,
    hospitalization_rate = .2,
    hospitalization_period = 7,
    prop_vaccinated = 1 - 1 / 15,
    quarantine_period = 21,
    quarantine_willingness = 1,
    isolation_period = 4,
    ...
    ) {
  # Check input parameters
  stopifnot_int(n, lb = 1)
  stopifnot_int(prevalence, lb = 0, ub = n)
  stopifnot_double(contact_rate, lb = 1e-10, ub = n)
  stopifnot_double(transmission_rate, lb = 0, ub = 1)
  stopifnot_double(vax_efficacy, lb = 0, ub = 1)
  # stopifnot_double(vax_improved_recovery, lb = 0, ub = 1)
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

  if ("vax_improved_recovery" %in% names(list(...)))
    warning(
      "The argument 'vax_improved_recovery' is no longer used and has been ",
      "removed. Future versions may not accept it.",
      call. = FALSE,
      immediate. = TRUE
    )

  structure(
    ModelMeaslesSchool_cpp(
      n,
      prevalence,
      contact_rate,
      transmission_rate,
      vax_efficacy,
      0.0, # vax_improved_recovery,
      incubation_period,
      prodromal_period,
      rash_period,
      days_undetected,
      hospitalization_rate,
      hospitalization_period,
      prop_vaccinated,
      quarantine_period,
      quarantine_willingness,
      isolation_period
    ),
    class = c("epiworld_measlesquarantine", "epiworld_model")
  )

}

#' @export
#' @rdname epiworldR-deprecated
ModelMeaslesQuarantine <- function(...) {
  .Deprecated(
    "ModelMeaslesSchool",
    package = "epiworldR",
    msg =
      "ModelMeaslesQuarantine is deprecated. Use ModelMeaslesSchool instead."
  )

  ModelMeaslesSchool(...)

}
