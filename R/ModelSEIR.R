#' Susceptible Exposed Infected Recovered model (SEIR)
#'
#' @param name String. Name of the virus.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus. 
#' @param x Object of class SEIR. 
#' @param ... Currently ignore.
#' @export
#' @family Models
#' @aliases epiworld_seir
#' @returns
#' - The `ModelSEIR`function returns a model of class [epiworld_model].
#' @examples 
#' model_seir <- ModelSEIR(name = "COVID-19", prevalence = 0.01, 
#' transmission_rate = 0.9, recovery_rate = 0.1, incubation_days = 4)
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
#' 
#' plot(model_seir, main = "SEIR Model")
#' @seealso epiworld-methods
ModelSEIR <- function(
    name, prevalence, transmission_rate, incubation_days, recovery_rate
) {
  
  structure(
    ModelSEIR_cpp(name, prevalence, transmission_rate, incubation_days, recovery_rate),
    class = c("epiworld_seir", "epiworld_model")
  )
  
}

#' @rdname ModelSEIR
#' @param main Title of the plot
#' @returns The `plot` function returns a plot of the SEIR model of class 
#' [epiworld_model].
#' @export
plot.epiworld_seir <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
