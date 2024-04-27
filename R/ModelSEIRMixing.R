#' Susceptible Exposed Infected Removed model (SEIR) with mixing
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of 
#' transmission.
#' @param incubation_days Numeric scalar. Average number of days in the
#' incubation period.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery.
#' @param x Object of class SIRCONN. 
#' @param ... Currently ignore.
#' @param n Number of individuals in the population.
#' @param contact_matrix Matrix of contact rates between individuals.
#' @export
#' @family Models
#' @details 
#' The `contact_matrix` is a matrix of contact rates between entities. The
#' matrix should be of size `n x n`, where `n` is the number of entities. 
#' This is a row-stochastic matrix, i.e., the sum of each row should be 1.
#' 
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @returns
#' - The `ModelSEIRMixing`function returns a model of class [epiworld_model].
#' @aliases epiworld_seirmixing
#' 
#' @examples 
#' # TBD
#' @seealso epiworld-methods
ModelSEIRMixing <- function(
    name, n, prevalence, contact_rate, transmission_rate, 
    incubation_days, recovery_rate, contact_matrix
) {
  
  structure(
    ModelSEIRCONN_cpp(
      name, n, prevalence, contact_rate, 
      transmission_rate, incubation_days, 
      recovery_rate, as.vector(contact_matrix)
      ),
    class = c("epiworld_seirmixing", "epiworld_model")
  )
  
}

#' @rdname ModelSEIRMixing
#' @export
#' @returns The `plot` function returns a plot of the SEIRMixing model of class 
#' [epiworld_model].
#' @param main Title of the plot
plot.epiworld_seirmixing <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
}
