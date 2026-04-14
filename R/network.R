#' Network models in epiworld
#' @name network-models
NULL

#' @param model Model object of class [epiworld_model].
#' @param block_sizes Integer vector describing the size of each block in the SBM.
#' @param mixing_matrix Numeric matrix describing the mixing between blocks in the SBM.
#' @return
#' - The function `agents_sbm` generates an stochastic block model (SBM) network.
#' @rdname network-models
#' @export
agents_sbm <- function(model, block_sizes, mixing_matrix) UseMethod("agents_sbm")

#' @rdname network-models
#' @examples
#' # Initializing SIR model with SBM network
#' sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, transmission_rate = 0.9,
#'   recovery_rate = 0.1)
#' agents_sbm(
#'   sir,
#'   block_sizes = c(500, 500),
#'   mixing_matrix = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
#' )
#' run(sir, ndays = 100, seed = 1912)
#' sir
#' @export
agents_sbm.epiworld_model <- function(model, block_sizes, mixing_matrix) {
  agents_sbm_cpp(
    model,
    as.integer(block_sizes),
    as.vector(mixing_matrix),
    row_major = FALSE
  )
  invisible(model)
}
