#' Susceptible Infected Removed model (SIR connected)
#' @param name Name of the virus
#' @param prevalence Initial proportion of individuals with the virus.
#' @param reproductive_number Numeric scalar. Reproductive number.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of transmission.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param n,m,days,seed,x,... to be documented
#' @export
#' @family Models
#' @aliases epiworld_sirconn
#' @examples 
#' model_sirconn <- ModelSIRCONN(name = "COVID-19", prevalence = 0.01, reproductive_number = 5, prob_transmission = 0.4, prob_recovery = 0.95)
#' # Adding a small world population
#' agents_smallworld(
#'   model_sirconn,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Initializing
#' init(model_sirconn, days = 100, seed = 1912)
#' # Running and printing
#' run(model_sirconn)
#' model_sirconn
#' 
#' 
#' 
#' 
#' 
#' 
#' 
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
 plot_epi(x, main = "SIRCONN Model", counts_scale, ...)
}