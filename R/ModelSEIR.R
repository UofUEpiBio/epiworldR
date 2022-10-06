#' Susceptible Infected Susceptible model (SEIR)
#'
#' @param name Name of the virus
#' @param prevalence a number
#' @param infectiousness a number
#' @param incubation_days a number
#' @param recovery a number
#' @param x to be documented
#' @param ... to be documented
#' @export
#' @family Models
#' @aliases epiworld_seir
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
