#' SIRD model
#' @param name String. Name of the virus.
#'
#' Susceptible-Infected-Recovered-Deceased model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1166587407#The_SIRD_model)).
#'
#' @param name String. Name of the virus

#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @param death_rate Numeric scalar between 0 and 1. Rate of death from virus.
#' @export
#' @concept general-models
#' @details
#' [set_initial()] places agents, at the start of each run, in stages without
#' the virus (e.g., `set_initial(model, assign_shares(c(Recovered = 0.2)))`); the
#' model's stages are listed by [get_stages()]. The model is built by
#' [make_model()] (see [get_model_info()] for its arguments' limits).
#' @section Model diagram:
#' ![](sird.png "SIRD Diagram")
#' @family Models
#' @aliases epiworld_sird
#' @returns
#' - The `ModelSIRD` function returns a model of class [epiworld_model].
#' @examples
#' model_sird <- ModelSIRD(
#'   name              = "COVID-19",
#'   n                 = 1000,
#'   prevalence        = 0.01,
#'   transmission_prob = 0.9,
#'   recovery_rate     = 0.1,
#'   death_rate        = 0.01
#' )
#'
#' # Adding a small world population
#' add_layer(model_sird, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_sird, ndays = 100, seed = 1912)
#' model_sird
#'
#' # Plotting
#' plot(model_sird)
#' @seealso epiworld-methods
ModelSIRD <- function(
  name, n, prevalence, transmission_prob, recovery_rate, death_rate
) {
  # The registry checks the values (see get_model_info())

  make_model(
    "SIRD",
    list(
      n                 = n,
      prevalence        = prevalence,
      transmission_prob = transmission_prob,
      recovery_rate     = recovery_rate,
      death_rate        = death_rate
    ),
    virus_name = name
  )

}
