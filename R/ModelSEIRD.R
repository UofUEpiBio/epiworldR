#' Susceptible-Exposed-Infected-Recovered-Deceased model (SEIRD)
#'
#' @param name String. Name of the virus.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_rate Numeric scalar between 0 and 1. Virus's rate of 
#' infection.
#' @param incubation_days Numeric scalar greater than 0. Average number of 
#' incubation days.
#' @param recovery_rate Numeric scalar between 0 and 1. Rate of recovery_rate from virus. 
#' @param death_rate Numeric scalar between 0 and 1. Rate of death from virus.
#' @param x Object of class SEIRD. 
#' @param ... Currently ignore.
#' @export
#' @details 
#' The [initial_state] function allows the user to set the initial state of the
#' model. The user must provide a vector of proportions indicating the following
#' values: (1) Proportion of exposed agents who are infected, (2)
#' proportion of non-infected agents already removed, and (3) proportion of
#' non-ifected agents already deceased.
#' @family Models
#' @aliases epiworld_seird
#' @returns
#' - The `ModelSEIRD`function returns a model of class [epiworld_model].
#' @examples 
#' model_seird <- ModelSEIRD(name = "COVID-19", prevalence = 0.01, 
#' transmission_rate = 0.9, recovery_rate = 0.1, incubation_days = 4, 
#' death_rate = 0.01)
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   model_seird,
#'   n = 100000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Running and printing
#' run(model_seird, ndays = 100, seed = 1912)
#' model_seird
#' 
#' plot(model_seird, main = "SEIRD Model")
#' @seealso epiworld-methods
ModelSEIRD <- function(
    name,
    prevalence,
    transmission_rate,
    incubation_days,
    recovery_rate,
    death_rate
) {
  
  structure(
    ModelSEIRD_cpp(
      name              = name,
      prevalence        = prevalence,
      transmission_rate = transmission_rate,
      incubation_days   = incubation_days,
      recovery_rate     = recovery_rate,
      death_rate        = death_rate
      ),
    class = c("epiworld_seird", "epiworld_model")
  )
  
}

#' @rdname ModelSEIRD
#' @param main Title of the plot
#' @returns The `plot` function returns a plot of the SEIRD model of class 
#' [epiworld_model].
#' @export
plot.epiworld_seird <- function(x, main = get_name(x), ...) { # col = NULL
 plot_epi(x, main = main, ...)
} 
