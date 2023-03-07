#' Susceptible Infected Susceptible model (SEIR)
#'
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param infectiousness Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery Numeric scalar between 0 and 1. Rate of recovery from virus. 
#' @param x Object of class SEIR. 
#' @param ... Currently ignore.
#' @export
#' @family Models
#' @aliases epiworld_seir
#' @examples 
#' model_seir <- ModelSEIR(name = "COVID-19", prevalence = 0.01, 
#' infectiousness = 0.9, recovery = 0.1, incubation_days = 4)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_seir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Running and printing
#' run(model_seir, ndays = 100, seed = 1912)
#' model_seir
#' @seealso epiworld-methods
ModelSEIR <- function(
    name, prevalence, infectiousness, incubation_days, recovery
) {
  
  structure(
    ModelSEIR_cpp(name, prevalence, infectiousness, incubation_days, recovery),
    class = c("epiworld_model", "epiworld_seir")
  )
  
}

#' @rdname ModelSEIR
#' @param main Title of the plot
#' @export
plot.epiworld_seir <- function(x, main = "SEIR Model", ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
