#' Susceptible-Exposed-Infected-Recovered-Deceased model (SEIRD)
#'
#' @param name String. Name of the virus.
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param latent_period Numeric scalar greater than 0. Average number of
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @param death_rate Numeric scalar between 0 and 1. Rate of death from virus.
#' @export
#' @concept general-models
#' @details
#' [set_initial()] places agents, at the start of each run, in stages without
#' the virus (e.g., `set_initial(model, assign_shares(c(Removed = 0.2)))`); the
#' model's stages are listed by [get_stages()]. The model is built by
#' [make_model()] (see [get_model_info()] for its arguments' limits).
#' @section Model diagram:
#' ![](seird.png "SEIRD Diagram")
#' @family Models
#' @aliases epiworld_seird
#' @returns
#' - The `ModelSEIRD`function returns a model of class [epiworld_model].
#' @examples
#' model_seird <- ModelSEIRD(name = "COVID-19", n = 100000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1, latent_period = 4,
#'   death_rate = 0.01)
#'
#' # Adding a small world population
#' add_layer(model_seird, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_seird, ndays = 100, seed = 1912)
#' model_seird
#'
#' plot(model_seird, main = "SEIRD Model")
#' @seealso epiworld-methods
ModelSEIRD <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  latent_period,
  recovery_rate,
  death_rate
) {
  # The registry checks the values (see get_model_info())

  make_model(
    "SEIRD",
    list(
      n                 = n,
      prevalence        = prevalence,
      transmission_prob = transmission_prob,
      latent_period     = latent_period,
      recovery_rate     = recovery_rate,
      death_rate        = death_rate
    ),
    virus_name = name
  )

}
