#' Susceptible Exposed Infected Removed Deceased model (SEIRD connected)
#' 
#' The SEIRD connected model implements a model where all agents are connected.
#' This is equivalent to a compartmental model ([wiki](https://en.wikipedia.org/w/index.php?title=Compartmental_models_in_epidemiology&oldid=1155757336#The_SEIR_model)).
#' 
#' @param name String. Name of the virus.
#' @param n Integer greater than zero. Population size.
#' @param prevalence Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery_rate.
#' @param death_rate Numeric scalar between 0 and 1. Probability of death.
#' @param x Object of class SEIRCONN. 
#' @param ... Currently ignore. 
#' @param n Number of individuals in the population.
#' @export
#' @family Models
#' @aliases epiworld_seirconn
#' @returns
#' - The `ModelSEIRDCONN`function returns a model of class [epiworld_model].
#' @examples 
#' # An example with COVID-19
#' model_seirdconn <- ModelSEIRDCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01, 
#'   n                   = 10000,
#'   contact_rate        = 2, 
#'   incubation_days     = 7, 
#'   transmission_rate   = 0.5,
#'   recovery_rate       = 0.3,
#'   death_rate          = 0.01
#' )
#'   
#' # Running and printing
#' run(model_seirdconn, ndays = 100, seed = 1912)
#' model_seirdconn
#' 
#' plot(model_seirdconn)
#' 
#' # Adding the flu
#' flu <- virus("Flu", .9, 1/7)
#' add_virus(model_seirdconn, flu, .001)
#' 
#' #' # Running and printing
#' run(model_seirdconn, ndays = 100, seed = 1912)
#' model_seirdconn
#' 
#' plot(model_seirdconn)
#' @seealso epiworld-methods
ModelSEIRDCONN <- function(
    name, n, prevalence, contact_rate, transmission_rate, 
    incubation_days, recovery_rate, death_rate
) {
  
  structure(
    ModelSEIRDCONN_cpp(name, n, prevalence, contact_rate, 
                      transmission_rate, incubation_days, recovery_rate,
                      death_rate),
    class = c("epiworld_seirdconn", "epiworld_model")
  )
  
}

#' @rdname ModelSEIRDCONN
#' @export
#' @returns The `plot` function returns a plot of the SEIRDCONN model of class 
#' [epiworld_model].
#' @param main Title of the plot.
plot.epiworld_seirdconn <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
