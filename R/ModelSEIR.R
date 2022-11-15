#' Susceptible Infected Susceptible model (SEIR)
#'
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param x to be documented
#' @param ... to be documented
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

#' @param m to be documented
#'
#' @param days to be documented
#' @param seed to be documented
#'
#' @rdname ModelSEIR
#' @export
init.epiworld_seir <- function(m, days, seed) {
  init_sir(m, days, seed)
}

#' @rdname ModelSEIR
#' @export
print.epiworld_seir <- function(x, ...) {
  print_sir(x)
}

#' @param n to be documented
#' @param k to be documented
#' @param d to be documented
#' @param p to be documented
#'
#' @rdname ModelSEIR
#' @export
agents_smallworld.epiworld_seir <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSEIR
#' @export
run.epiworld_seir <- function(m) {
  run_sir(m)
}

#' @rdname ModelSEIR
#' @export
plot.epiworld_seir <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SEIR Model", counts_scale, ...)
}
