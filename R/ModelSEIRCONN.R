#' Susceptible Exposed Infected Removed model (SEIR connected)
#' 
#' The SEIR connected model implements a model where all agents are connected.
#' This is equivalent to a compartmental model.
#' 
#' @param name String. Name of the virus
#' @param n Integer greater than zero. Population size.
#' @param prevalence Initial proportion of individuals with the virus.
#' @param reproductive_number Numeric scalar. Reproductive number.
#' @param prob_transmission Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param prob_recovery Numeric scalar between 0 and 1. Probability of recovery.
#' @param m  
#' @param days Numeric scalar. Number of days the simulation is to run for. 
#' @param seed Seed to set for initializing random number generator.
#' @param x 
#' @param ...
#' @param n Number of individuals in the population.
#' @param k
#' @param d 
#' @param p 
#' @export
#' @family Models
#' @aliases epiworld_seirconn
#' @examples 
#' model_seirconn <- ModelSEIRCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01, 
#'   n                   = 10000,
#'   reproductive_number = 4, 
#'   incubation_days     = 7, 
#'   prob_transmission   = 0.5,
#'   prob_recovery       = 0.99
#' )
#'   
#' # Initializing
#' init(model_seirconn, days = 100, seed = 1912)
#' # Running and printing
#' run(model_seirconn)
#' model_seirconn

ModelSEIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, 
    incubation_days, prob_recovery
) {
  
  structure(
    ModelSEIRCONN_cpp(name, n, prevalence, reproductive_number, 
                      prob_transmission, incubation_days, prob_recovery),
    class = "epiworld_seirconn"
  )
  
}

#' @rdname ModelSEIRCONN
#' @export
init.epiworld_seirconn <- function(m, days, seed) {
  init_cpp(m, days, seed)
}

#' @rdname ModelSEIRCONN
#' @export
print.epiworld_seirconn <- function(x, ...) {
  print_cpp(x)
}

#' #' @rdname ModelSEIRCONN
#' #' @export
#' agents_smallworld.epiworld_seirconn <- function(m, n, k, d, p) {
#'   agents_smallworld_cpp(m, n, k, d, p)
#' }

#' @rdname ModelSEIRCONN
#' @export
run.epiworld_seirconn <- function(m) {
  run_cpp(m)
}

#' @rdname ModelSEIRCONN
#' @export
plot.epiworld_seirconn <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SEIRCONN Model", ...)
}
