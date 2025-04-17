#' Susceptible Infected Removed model (SIR) with mixing
#' @param name String. Name of the virus
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param contact_rate Numeric scalar. Average number of contacts per step.
#' @param transmission_rate Numeric scalar between 0 and 1. Probability of
#' transmission.
#' @param recovery_rate Numeric scalar between 0 and 1. Probability of recovery.
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
#' - The `ModelSIRMixing`function returns a model of class [epiworld_model].
#' @aliases epiworld_sirmixing
#'
#' @examples
#' # From the vignette
#'
#' # Start off creating three entities.
#' # Individuals will be distribured randomly between the three.
#' e1 <- entity("Population 1", 3e3, as_proportion = FALSE)
#' e2 <- entity("Population 2", 3e3, as_proportion = FALSE)
#' e3 <- entity("Population 3", 3e3, as_proportion = FALSE)
#'
#' # Row-stochastic matrix (rowsums 1)
#' cmatrix <- c(
#'   c(0.9, 0.05, 0.05),
#'   c(0.1, 0.8, 0.1),
#'   c(0.1, 0.2, 0.7)
#' ) |> matrix(byrow = TRUE, nrow = 3)
#'
#' N <- 9e3
#'
#' flu_model <- ModelSIRMixing(
#'   name              = "Flu",
#'   n                 = N,
#'   prevalence        = 1 / N,
#'   contact_rate      = 20,
#'   transmission_rate = 0.1,
#'   recovery_rate     = 1 / 7,
#'   contact_matrix    = cmatrix
#' )
#'
#' # Adding the entities to the model
#' flu_model |>
#'   add_entity(e1) |>
#'   add_entity(e2) |>
#'   add_entity(e3)
#'
#' set.seed(331)
#' run(flu_model, ndays = 100)
#' summary(flu_model)
#' plot_incidence(flu_model)
#'
#' @seealso epiworld-methods
ModelSIRMixing <- function(
    name,
    n,
    prevalence,
    contact_rate,
    transmission_rate,
    recovery_rate,
    contact_matrix
    ) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n)
  stopifnot_double(prevalence)
  stopifnot_double(contact_rate)
  stopifnot_double(transmission_rate)
  stopifnot_double(recovery_rate)
  stopifany_na(contact_matrix)

  structure(
    ModelSIRMixing_cpp(
      name, n, prevalence, contact_rate,
      transmission_rate, recovery_rate, as.vector(contact_matrix)
    ),
    class = c("epiworld_sirmixing", "epiworld_model")
  )

}
