#' Susceptible Infected Susceptible model (SIS)
#' @param name Name of the virus
#' @param prevalence a number
#' @param infectiousness a number
#' @param recovery a number
#' @export
#' @family Models
ModelSIS <- function(
    name, prevalence, infectiousness, recovery
) {
  
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
