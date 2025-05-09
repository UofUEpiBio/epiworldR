#' SURV model
#'
#'
#' @param name String. Name of the virus.
#' @param prevalence Initial number of individuals with the virus.
#' @param efficacy_vax Double. Efficacy of the vaccine.
#' (1 - P(acquire the disease)).
#' @param latent_period Double. Shape parameter of a 'Gamma(latent_period, 1)'
#' distribution. This coincides with the expected number of latent days.
#' @param infect_period Double. Shape parameter of a 'Gamma(infected_period, 1)'
#' distribution. This coincides with the expected number of infectious days.
#' @param prob_symptoms Double. Probability of generating symptoms.
#' @param prop_vaccinated Double. Probability of vaccination. Coincides with
#' the initial prevalence of vaccinated individuals.
#' @param prop_vax_redux_transm Double. Factor by which the vaccine reduces
#' transmissibility.
#' @param prop_vax_redux_infect Double. Factor by which the vaccine reduces
#' the chances of becoming infected.
#' @param surveillance_prob Double. Probability of testing an agent.
#' @param transmission_rate Double. Raw transmission probability.
#' @param prob_death Double. Raw probability of death for symptomatic
#' individuals.
#' @param prob_noreinfect Double. Probability of no re-infection.
#' @export
#' @family Models
#' @aliases epiworld_surv
#' @returns
#' - The `ModelSURV`function returns a model of class [epiworld_model].
#' @examples
#' model_surv <- ModelSURV(
#'   name                  = "COVID-19",
#'   prevalence            = 20,
#'   efficacy_vax          = 0.6,
#'   latent_period         = 4,
#'   infect_period         = 5,
#'   prob_symptoms         = 0.5,
#'   prop_vaccinated       = 0.7,
#'   prop_vax_redux_transm = 0.8,
#'   prop_vax_redux_infect = 0.95,
#'   surveillance_prob     = 0.1,
#'   transmission_rate     = 0.2,
#'   prob_death            = 0.001,
#'   prob_noreinfect       = 0.5
#' )
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_surv,
#'   n = 10000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
#'
#' # Running and printing
#' run(model_surv, ndays = 100, seed = 1912)
#' model_surv
#'
#' # Plotting
#' plot(model_surv, main = "SURV Model")
#'
#' @seealso epiworld-methods
ModelSURV <- function(
    name,
    prevalence,
    efficacy_vax,
    latent_period,
    infect_period,
    prob_symptoms,
    prop_vaccinated,
    prop_vax_redux_transm,
    prop_vax_redux_infect,
    surveillance_prob,
    transmission_rate,
    prob_death,
    prob_noreinfect
    ) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(efficacy_vax)
  stopifnot_double(latent_period)
  stopifnot_double(infect_period)
  stopifnot_double(prob_symptoms)
  stopifnot_double(prop_vaccinated)
  stopifnot_double(prop_vax_redux_transm)
  stopifnot_double(prop_vax_redux_infect)
  stopifnot_double(surveillance_prob)
  stopifnot_double(transmission_rate)
  stopifnot_double(prob_death)
  stopifnot_double(prob_noreinfect)

  structure(
    ModelSURV_cpp(name, prevalence, efficacy_vax, latent_period, infect_period,
      prob_symptoms, prop_vaccinated, prop_vax_redux_transm,
      prop_vax_redux_infect, surveillance_prob, transmission_rate,
      prob_death, prob_noreinfect),
    class = c("epiworld_surv", "epiworld_model")
  )

}
