#' Susceptible Exposed Infected Removed model (SEIR connected)
#' 
#' The SEIR connected model implements a model where all agents are connected.
#' This is equivalent to a compartmental model.
#' 
#' @param name Name of the virus
#' @param n Integer greater than zero. Population size.
#' @param prevalence Initial proportion of individuals with the virus.
#' @param reproductive_number Numeric scalar. Reproductive number.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of transmission.
#' @param incubation_days Numeric scalar greater than 0. Average number of incubation days.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param m  to be documented
#' @param days to be documented
#' @param seed to be documented
#' @param x to be documented
#' @param ... to be documented
#' @export
#' @family Models
#' @aliases epiworld_seirconn
ModelSEIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, incubation_days, prob_recovery
) {
  
  structure(
    ModelSEIRCONN_cpp(name, n, prevalence, reproductive_number, prob_transmission, incubation_days, prob_recovery),
    class = "epiworld_seirconn"
  )
  
}

#' @rdname ModelSEIRCONN
#' @export
init.epiworld_seirconn <- function(m, days, seed) {
  init_seirconn(m, days, seed)
}

#' @rdname ModelSEIRCONN
#' @export
print.epiworld_seirconn <- function(x, ...) {
  print_seirconn(x)
}

#' #' @rdname ModelSEIRCONN
#' #' @export
#' agents_smallworld.epiworld_seirconn <- function(m, n, k, d, p) {
#'   agents_smallworld_sir(m, n, k, d, p)
#' }

#' @rdname ModelSEIRCONN
#' @export
run.epiworld_seirconn <- function(m) {
  run_seirconn(m)
}

#' @rdname ModelSEIRCONN
#' @export
plot.epiworld_seirconn <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SEIRCONN Model", counts_scale, ...)
}
