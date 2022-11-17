#' SIR model
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param x
#' @param ...
#' @param n Number of individuals in the population.
#' @param k
#' @param d 
#' @param p 
#' @export
#' @family Models
#' @aliases epiworld_sir
#' @examples 
#' model_sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, 
#'                       infectiousness = 0.9, recovery = 0.1)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_sir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Initializing
#' init(model_sir, days = 100, seed = 1912)
#' # Running and printing
#' run(model_sir)
#' model_sir

ModelSIR <- function(
    name, prevalence, infectiousness, recovery
) {
  
  structure(
    ModelSIR_cpp(name, prevalence, infectiousness, recovery),
    class = "epiworld_sir"
  )
  
}

#' @rdname ModelSIR
#' @export
init.epiworld_sir <- function(m, days, seed) {
  init_sir(m, days, seed)
}

#' @rdname ModelSIR
#' @export
print.epiworld_sir <- function(x, ...) {
  print_sir(x)
}

#' @rdname ModelSIR
#' @export
agents_smallworld.epiworld_sir <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSIR
#' @export
run.epiworld_sir <- function(m) {
  run_sir(m)
}

#' @rdname ModelSIR
#' @export
plot.epiworld_sir <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SIR Model", counts_scale, ...)
}

