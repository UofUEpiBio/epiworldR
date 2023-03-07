#' Susceptible Infected Removed model (SIR connected)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param reproductive_number Numeric scalar. Reproductive number.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param x Object of class SIRCONN. 
#' @param ... Currently ignore.
#' @param n Number of individuals in the population.
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
#' # Running and printing
#' run(model_sirconn, ndays = 100, seed = 1912)
#' model_sirconn
#' @seealso epiworld-methods
ModelSIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, prob_recovery
) {
  
  structure(
    ModelSIRCONN_cpp(name, n, prevalence, reproductive_number, 
                     prob_transmission, prob_recovery),
    class = c("epiworld_model", "epiworld_sirconn")
  )
  
}

#' @rdname ModelSIRCONN
#' @export
#' @param main Title of the plot
plot.epiworld_sirconn <- function(x, main = "SIRCONN Model", ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
