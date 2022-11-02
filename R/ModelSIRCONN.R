#' Susceptible Infected Removed model (SIR connected)
#' @param name Name of the virus
#' @param prevalence a number
#' @param reproductive_number a number
#' @param prob_transmission a number
#' @param prob_recovery a number
#' @param n,m,days,seed,x,... to be documented
#' @export
#' @family Models
#' @aliases epiworld_sirconn
ModelSIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, prob_recovery
) {
  
  structure(
    ModelSIRCONN_cpp(name, n, prevalence, reproductive_number, prob_transmission, prob_recovery),
    class = "epiworld_sirconn"
  )
  
}

#' @rdname ModelSIRCONN
#' @export
init.epiworld_sirconn <- function(m, days, seed) {
  init_sirconn(m, days, seed)
}

#' @rdname ModelSIRCONN
#' @export
print.epiworld_sirconn <- function(x, ...) {
  print_sirconn(x)
}

#' @rdname ModelSIRCONN
#' @export
run.epiworld_sirconn <- function(m) {
  run_sirconn(m)
}

#' @rdname ModelSIRCONN
#' @export
plot.epiworld_sirconn <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SIRCONN Model", ...)
}