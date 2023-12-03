#' Susceptible Infected Removed model (SIR connected)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery.
#' @param x Object of class SIRCONN. 
#' @param ... Currently ignore.
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @details 
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
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
#'   transmission_rate   = 0.4,
#'   recovery_rate       = 0.95
#' )
#'   
#' # Running and printing
#' run(model_sirconn, ndays = 100, seed = 1912)
#' model_sirconn
#' 
#' plot(model_sirconn,  main = "SIRCONN Model")
#' @seealso epiworld-methods
ModelSIRCONN <- function(
    name, n, prevalence, contact_rate, transmission_rate, recovery_rate
) {
  
  structure(
    ModelSIRCONN_cpp(name, n, prevalence, contact_rate, 
                     transmission_rate, recovery_rate),
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
