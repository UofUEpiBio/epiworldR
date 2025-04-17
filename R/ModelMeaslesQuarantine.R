#' Measles model with quarantine
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
#' @param hospitalization_duration Average number of days in hospital.
#' @param prop_vaccinated Proportion of the population vaccinated.
#' @param quarantine_days Total duration of quarantine.
#' @param quarantine_willingness Probability of accepting quarantine.
#' @details
#' The basic reproductive number in Measles is estimated to be about 15.
#' By default, the contact rate of the model is set so that the R0 matches
#' 15.
#' @export
#' @family Models
#' @aliases epiworld_measlesquarantine
#' @returns
#' - The `ModelSEIRCONN`function returns a model of class [epiworld_model].
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
ModelMeaslesQuarantine <- function(
    n,
    prevalence = 1,
    contact_rate = 15 / transmission_rate / (rash_period + prodromal_period),
    transmission_rate = .9,
    vax_efficacy = .99,
    vax_improved_recovery = .5,
    incubation_period = 12,
    prodromal_period = 4,
    rash_period = 3,
    days_undetected = 2,
    hospitalization_rate = .2,
    hospitalization_duration = 7,
    prop_vaccinated = 1 - 1 / 15,
    quarantine_days = 21,
    quarantine_willingness = 1
    ) {
  # Check input parameters
  stopifnot_int(n)
  stopifnot_int(prevalence)
  stopifnot_double(contact_rate)
  stopifnot_double(transmission_rate)
  stopifnot_double(vax_efficacy)
  stopifnot_double(vax_improved_recovery)
  stopifnot_double(incubation_period)
  stopifnot_double(prodromal_period)
  stopifnot_double(rash_period)
  stopifnot_double(days_undetected)
  stopifnot_double(hospitalization_rate)
  stopifnot_double(hospitalization_duration)
  stopifnot_double(prop_vaccinated)
  stopifnot_double(quarantine_days)
  stopifnot_double(quarantine_willingness)

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
      hospitalization_duration,
      prop_vaccinated,
      quarantine_days,
      quarantine_willingness
    ),
    class = c("epiworld_measlesquarantine", "epiworld_model")
  )

}

#' @rdname ModelMeaslesQuarantine
#' @export
#' @param x Object of class [epiworld_measlesquarantine].
#' @param main Title of the plot
#' @param ... Passed to [graphics::plot].
#' @returns The `plot` function returns a plot of the MeaslesQuarantine model of class
#' [epiworld_model].
plot.epiworld_measlesquarantine <- function(x, main = get_name(x), ...) { # col = NULL
  plot_epi(x, main = main, ...)
}
