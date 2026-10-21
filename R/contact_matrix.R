#' Get or set the contact matrix of a mixing layer
#' @param layer A mixing layer (class `epiworld_layer_mixing`), retrieved with
#' [get_layer()] (the mixing models name theirs `"mixing"`).
#' @param contact_matrix A numeric matrix with the contact matrix to be
#' set for the layer. The matrix should be square and have as many rows as the
#' layer has groups.
#'
#' @details
#' Set between runs, the contact matrix is kept for the following runs; set
#' during a run (e.g., from a global event), it lasts until the end of that
#' run.
#'
#' @return
#' - The `get_contact_matrix()` function returns the contact matrix of the
#' layer as a numeric matrix.
#' @examples
#' cmatrix <- (c(
#'   c(0.9, 0.05, 0.05),
#'   c(0.1, 0.8, 0.1),
#'   c(0.1, 0.2, 0.7)
#' ) * 20) |> matrix(byrow = TRUE, nrow = 3)
#'
#'
#' N <- 9e3
#'
#' flu_model <- ModelSIRMixing(
#'   name              = "Flu",
#'   n                 = N,
#'   prevalence        = 1 / N,
#'   transmission_rate = 0.1,
#'   recovery_rate     = 1 / 7,
#'   contact_matrix    = cmatrix
#' )
#'
#' mixing <- get_layer(flu_model, "mixing")
#' get_contact_matrix(mixing)
#'
#' # Modifying the contact matrix
#' new_cmatrix <- (c(
#'   c(0.8, 0.1, 0.1),
#'   c(0.1, 0.7, 0.2),
#'   c(0.1, 0.2, 0.7)
#' ) * 20) |> matrix(byrow = TRUE, nrow = 3)
#' set_contact_matrix(mixing, new_cmatrix)
#' get_contact_matrix(mixing)
#' @export
#' @concept mixing-models
get_contact_matrix <- function(layer) {
  stopifnot_layer(layer, "mixing")
  ans <- get_contact_matrix_cpp(layer)
  matrix(
    ans,
    nrow = sqrt(length(ans)),
    byrow = FALSE
  )

}

#' @rdname get_contact_matrix
#' @export
#' @return
#' - The `set_contact_matrix()` function sets the contact matrix of the layer
#' and returns the layer invisibly.
set_contact_matrix <- function(layer, contact_matrix) {
  stopifnot_layer(layer, "mixing")
  stopifany_na(contact_matrix)
  set_contact_matrix_cpp(layer, as.vector(contact_matrix))
  invisible(layer)
}
