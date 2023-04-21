#' Susceptible Infected Susceptible model (SIS)
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection 
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param x Object of class SIS. 
#' @param ... Currently ignore. 
#' @export
#' @family Models
#' @aliases epiworld_sis
#' @examples 
#' model_sis <- ModelSIS(name = "COVID-19", prevalence = 0.01, 
#'                      infectiousness = 0.9, recovery = 0.1)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_sis,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Running and printing
#' run(model_sis, ndays = 100, seed = 1912)
#' model_sis
#' @seealso epiworld-methods
ModelSIS <- function(
    name, prevalence, infectiousness, recovery) {
  
  structure(
    ModelSIS_cpp(name, prevalence, infectiousness, recovery),
    class = c("epiworld_sis", "epiworld_model")
  )
  
}


#' @rdname ModelSIS
#' @export
#' @param main Title of the plot
plot.epiworld_sis <- function(x, main = "SIS Model",...) { # col = NULL
 plot_epi(x, main = main, ...)
}
