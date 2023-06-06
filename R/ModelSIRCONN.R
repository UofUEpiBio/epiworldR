#' Susceptible Infected Removed model (SIR connected)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param x Object of class SIRCONN. 
#' @param ... Currently ignore.
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @returns
#' - The `ModelSIRCONN`function returns a model of class [epiworld_model].
#' @aliases epiworld_sirconn
#' 
#' @examples 
#' model_sirconn <- ModelSIRCONN(
#'   name                = "COVID-19",
#'   n                   = 10000,
#'   prevalence          = 0.01,
#'   contact_rate        = 5,
#'   prob_transmission   = 0.4,
#'   prob_recovery       = 0.95
#' )
#'   
#' # Running and printing
#' run(model_sirconn, ndays = 100, seed = 1912)
#' model_sirconn
#' 
#' plot(model_sirconn,  main = "SIRCONN Model")
#' @seealso epiworld-methods
ModelSIRCONN <- function(
    name, n, prevalence, contact_rate, prob_transmission, prob_recovery
) {
  
  structure(
    ModelSIRCONN_cpp(name, n, prevalence, contact_rate, 
                     prob_transmission, prob_recovery),
    class = c("epiworld_sirconn", "epiworld_model")
  )
  
}

#' @rdname ModelSIRCONN
#' @export
#' @returns The `plot` function returns a plot of the SIRCONN model of class 
#' [epiworld_model].
#' @param main Title of the plot
plot.epiworld_sirconn <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
