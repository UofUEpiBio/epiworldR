#' Susceptible Infected Removed Deceased model (SIRD connected)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery.
#' @param death_rate Numeric scalar between 0 and 1. Probability of death.
#' @param x Object of class SIRDCONN. 
#' @param ... Currently ignore.
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @returns
#' - The `ModelSIRDCONN`function returns a model of class [epiworld_model].
#' @aliases epiworld_sirdconn
#' 
#' @examples 
#' model_sirdconn <- ModelSIRDCONN(
#'   name                = "COVID-19",
#'   n                   = 100000,
#'   prevalence          = 0.01,
#'   contact_rate        = 5,
#'   transmission_rate   = 0.4,
#'   recovery_rate       = 0.5,
#'   death_rate          = 0.1
#' )
#'   
#' # Running and printing
#' run(model_sirdconn, ndays = 100, seed = 1912)
#' model_sirdconn
#' 
#' plot(model_sirdconn,  main = "SIRDCONN Model")
#' @seealso epiworld-methods
ModelSIRDCONN <- function(
    name, n, prevalence, contact_rate, transmission_rate, recovery_rate,
    death_rate
) {
  
  structure(
    ModelSIRDCONN_cpp(name, n, prevalence, contact_rate, 
                     transmission_rate, recovery_rate, death_rate),
    class = c("epiworld_sirdconn", "epiworld_model")
  )
  
}

#' @rdname ModelSIRDCONN
#' @export
#' @returns The `plot` function returns a plot of the SIRDCONN model of class 
#' [epiworld_model].
#' @param main Title of the plot
plot.epiworld_sirdconn <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
