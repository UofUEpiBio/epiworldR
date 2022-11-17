#' Susceptible Infected Removed model (SIR connected)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param reproductive_number Numeric scalar. Reproductive number.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param n,m,days,seed,x,... to be documented
#' @export
#' @family Models
#' @aliases epiworld_sirconn
#' @examples 
#' model_sirconn <- ModelSIRCONN(
#'   name                = "COVID-19",
#'   n                   = 10000,
#'   prevalence          = 0.01,
#'   reproductive_number = 5,
#'   prob_transmission   = 0.4,
#'   prob_recovery       = 0.95
#' )
#'   
#' # Initializing
#' init(model_sirconn, days = 100, seed = 1912)
#' # Running and printing
#' run(model_sirconn)
#' model_sirconn
#' 
ModelSIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, prob_recovery
) {
  
  structure(
    ModelSIRCONN_cpp(name, n, prevalence, reproductive_number, 
                     prob_transmission, prob_recovery),
    class = "epiworld_sirconn"
  )
  
}

#' @rdname ModelSIRCONN
#' @export
init.epiworld_sirconn <- function(m, days, seed) {
  init_cpp(m, days, seed)
}

#' @rdname ModelSIRCONN
#' @export
print.epiworld_sirconn <- function(x, ...) {
  print_cpp(x)
}

#' @rdname ModelSIRCONN
#' @export
run.epiworld_sirconn <- function(m) {
  run_cpp(m)
}

#' @rdname ModelSIRCONN
#' @export
plot.epiworld_sirconn <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SIRCONN Model", ...)
}
