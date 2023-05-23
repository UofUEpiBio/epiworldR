#' SIR model
#' @param name String. Name of the virus.
#' 
#' Susceptible-Infected-Recovered model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SIR_model)).
#' 
#' @param name String. Name of the virus

#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param x Object of class SIR. 
#' @param ... Currently ignore. 
#' @export
#' @family Models
#' @aliases epiworld_sir
#' @examples 
#' model_sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, 
#'                       infectiousness = 0.9, recovery = 0.1)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_sir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Running and printing
#' run(model_sir, ndays = 100, seed = 1912)
#' model_sir
#' 
#' # Plotting
#' plot(model_sir)
#' @seealso epiworld-methods
ModelSIR <- function(
    name, prevalence, infectiousness, recovery
) {
  
  structure(
    ModelSIR_cpp(name, prevalence, infectiousness, recovery),
    class = c("epiworld_sir", "epiworld_model")
  )
  
}

#' @rdname ModelSIR
#' @export
#' @param main Title of the plot
plot.epiworld_sir <- function(x, main = "SIR Model", ...) { # col = NULL
 plot_epi(x, main = main, ...)
}

