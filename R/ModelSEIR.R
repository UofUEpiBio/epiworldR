#' Susceptible Exposed Infected Recovered model (SEIR)
#'
#' @param name String. Name of the virus.
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Virus's rate of
#' infection.
#' @param latent_period Numeric scalar greater than 0. Average number of
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus.
#' @export
#' @family Models
#' @concept general-models
#' @aliases epiworld_seir
#' @details
#' [set_initial()] places agents, at the start of each run, in stages without
#' the virus (e.g., `set_initial(model, assign_shares(c(Removed = 0.2)))`); the
#' model's stages are listed by [get_stages()]. The model is built by
#' [make_model()] (see [get_model_info()] for its arguments' limits).
#' @section Model diagram:
#' ![](seir.png "SEIR Diagram")
#' @returns
#' - The `ModelSEIR`function returns a model of class [epiworld_model].
#' @examples
#' model_seir <- ModelSEIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1, latent_period = 4)
#'
#' # Adding a small world population
#' add_layer(model_seir, layer_network("network", graph_smallworld(k = 5, p = .01)))
#'
#' # Running and printing
#' run(model_seir, ndays = 100, seed = 1912)
#' model_seir
#'
#' plot(model_seir, main = "SEIR Model")
#' @seealso epiworld-methods
ModelSEIR <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  latent_period,
  recovery_rate
) {
  # The registry checks the values (see get_model_info())

  make_model(
    "SEIR",
    list(
      n                 = n,
      prevalence        = prevalence,
      transmission_prob = transmission_prob,
      latent_period     = latent_period,
      recovery_rate     = recovery_rate
    ),
    virus_name = name
  )

}
