#' SIR model
#' @param name Name of the virus
#' @param prevalence a number
#' @param infectiousness a number
#' @param recovery a number
#' @param m,days,seed,x,...,n,k,d,p to be documented
#' @export
#' @family Models
#' @aliases epiworld_sir
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
