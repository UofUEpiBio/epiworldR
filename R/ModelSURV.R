#' Susceptible Infected Susceptible model (SURV)
#' @param name Name of the virus
#' @param prevalence a number
#' @param infectiousness a number
#' @param recovery a number
#' @export
#' @family Models
ModelSURV <- function(
    name, prevalence, infectiousness, recovery
) {
  
  structure(
    ModelSURV_cpp(name, prevalence, infectiousness, recovery),
    class = "epiworld_surv"
  )
  
}

#' @rdname ModelSURV
#' @export
init.epiworld_surv <- function(m, days, seed) {
  init_surv(m, days, seed)
}

#' @rdname ModelSURV
#' @export
print.epiworld_surv <- function(x, ...) {
  print_sir(x)
}

#' @rdname ModelSURV
#' @export
agents_smallworld.epiworld_surv <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSURV
#' @export
run.epiworld_surv <- function(m) {
  run_sir(m)
}
