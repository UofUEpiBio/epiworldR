#' Susceptible Exposed Infected Removed model (SEIR connected)
#' 
#' The SEIR connected model implements a model where all agents are connected.
#' This is equivalent to a compartmental model.
#' 
#' @param name String. Name of the virus
#' @param n Integer greater than zero. Population size.
#' @param prevalence Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param x Object of class SEIRCONN. 
#' @param ... Currently ignore. 
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @aliases epiworld_seirconn
#' @examples 
#' model_seirconn <- ModelSEIRCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01, 
#'   n                   = 10000,
#'   contact_rate = 4, 
#'   incubation_days     = 7, 
#'   prob_transmission   = 0.5,
#'   prob_recovery       = 0.99
#' )
#'   
#' # Running and printing
#' run(model_seirconn, ndays = 100, seed = 1912)
#' model_seirconn
#' @seealso epiworld-methods
ModelSEIRCONN <- function(
    name, n, prevalence, contact_rate, prob_transmission, 
    incubation_days, prob_recovery
) {
  
  structure(
    ModelSEIRCONN_cpp(name, n, prevalence, contact_rate, 
                      prob_transmission, incubation_days, prob_recovery),
    class = c("epiworld_seirconn", "epiworld_model")
  )
  
}

#' @rdname ModelSEIRCONN
#' @export
#' @param main Title of the plot
plot.epiworld_seirconn <- function(x, main = "SEIRCONN Model", ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
