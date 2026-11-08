#' SIR model
#' @param name String. Name of the virus.
#'
#' Susceptible-Infected-Recovered model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIR_model)).
#'
#' @param name String. Name of the virus

#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @export
#' @family Models
#' @concept general-models
#' @details
#' [set_initial()] places agents, at the start of each run, in stages without
#' the virus (e.g., `set_initial(model, assign_shares(c(Recovered = 0.2)))`); the
#' model's stages are listed by [get_stages()]. The model is built by
#' [make_model()] (see [get_model_info()] for its arguments' limits).
#' @section Model diagram:
#' ![](sir.png "SIR Diagram")
#' @aliases epiworld_sir
#' @returns
#' - The `ModelSIR` function returns a model of class [epiworld_model].
#' @examples
#' model_sir <- ModelSIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1)
#'
#' # Adding a small world population
#' add_layer(model_sir, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_sir, ndays = 100, seed = 1912)
#' model_sir
#'
#' # Plotting
#' plot(model_sir)
#' @seealso epiworld-methods
ModelSIR <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  recovery_rate
) {
  # The registry checks the values (see get_model_info())

  make_model(
    "SIR",
    list(
      n                 = n,
      prevalence        = prevalence,
      transmission_prob = transmission_prob,
      recovery_rate     = recovery_rate
    ),
    virus_name = name
  )

}
