#' Susceptible Infected Susceptible model (SIS)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection 
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param m,days,seed,x,...,n,k,d,p to be documented
#' @export
#' @family Models
#' @aliases epiworld_sis
#' @examples 
#' model_sis <- ModelSIS(name = "COVID-19", prevalence = 0.01, 
#'                      infectiousness = 0.9, recovery = 0.1)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_sis,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Initializing
#' init(model_sis, days = 100, seed = 1912)
#' # Running and printing
#' run(model_sis)
#' model_sis

ModelSIS <- function(
    name, prevalence, infectiousness, recovery) {
  
  structure(
    ModelSIS_cpp(name, prevalence, infectiousness, recovery),
    class = "epiworld_sis"
  )
  
}

#' @rdname ModelSIS
#' @export
init.epiworld_sis <- function(m, days, seed) {
  init_sis(m, days, seed)
}

#' @rdname ModelSIS
#' @export
print.epiworld_sis <- function(x, ...) {
  print_sir(x)
}

#' @rdname ModelSIS
#' @export
agents_smallworld.epiworld_sis <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSIS
#' @export
run.epiworld_sis <- function(m) {
  run_sir(m)
}

#' @rdname ModelSIS
#' @export
plot.epiworld_sis <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SIS Model", ...)
}
