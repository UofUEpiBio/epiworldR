#' Susceptible Infected Susceptible model (SEIR)
#'
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param m Model object.
#' @param days Number of days.
#' @param seed Seed to set for initializing random number generator.
#' @param x Object of class SEIR. 
#' @param ... Currently ignore. 
#' @param n Number of individuals in the population.
#' @param k Number of ties in the small world network.
#' @param d Logical scalar. Whether the graph is directed or not.
#' @param p Probability of rewiring.
#' @export
#' @family Models
#' @aliases epiworld_seir
#' @examples 
#' model_seir <- ModelSEIR(name = "COVID-19", prevalence = 0.01, 
#' infectiousness = 0.9, recovery = 0.1, incubation_days = 4)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_seir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Initializing
#' init(model_seir, days = 100, seed = 1912)
#' # Running and printing
#' run(model_seir)
#' model_seir

ModelSEIR <- function(
    name, prevalence, infectiousness, incubation_days, recovery
) {
  
  structure(
    ModelSEIR_cpp(name, prevalence, infectiousness, incubation_days, recovery),
    class = "epiworld_seir"
  )
  
}


#' @rdname ModelSEIR
#' @export
init.epiworld_seir <- function(m, days, seed) {
  init_cpp(m, days, seed)
}

#' @rdname ModelSEIR
#' @export
print.epiworld_seir <- function(x, ...) {
  print_cpp(x)
}

#' @rdname ModelSEIR
#' @export
agents_smallworld.epiworld_seir <- function(m, n, k, d, p) {
  agents_smallworld_cpp(m, n, k, d, p)
}

#' @rdname ModelSEIR
#' @export
run.epiworld_seir <- function(m) {
  run_cpp(m)
}

#' @rdname ModelSEIR
#' @export
plot.epiworld_seir <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SEIR Model", counts_scale, ...)
}
