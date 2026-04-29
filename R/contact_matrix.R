
#' Get or set the contact matrix of a mixing model
#' @param model An object of class `epiworld_model` that supports contact
#' matrices.
#' @param contact_matrix A numeric matrix with the contact matrix to be
#' set for the model. The matrix should be square and have the same number
#' of rows as the number of entities in the model.
#'
#' @return
#' - The `get_contact_matrix()` function returns the contact matrix of the
#' model as a numeric matrix.
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
#' get_contact_matrix(flu_model)
#' #' # Modifying the contact matrix
#' new_cmatrix <- (c(
#'   c(0.8, 0.1, 0.1),
#'   c(0.1, 0.7, 0.2),
#'   c(0.1, 0.2, 0.7)
#' ) * 20) |> matrix(byrow = TRUE, nrow = 3)
#' set_contact_matrix(flu_model, new_cmatrix)
#' get_contact_matrix(flu_model)
#' @export
#' @concept mixing-models
get_contact_matrix <- function(model) {
  stopifnot_model(model)
  ans <- get_contact_matrix_cpp(model)
  matrix(
    ans,
    nrow = sqrt(length(ans)),
    byrow = FALSE
  )

}

#' @rdname get_contact_matrix
#' @export
#' @param as_backup A logical value indicating whether to save the new
#' contact matrix as a backup in the model. If `TRUE` (default), the new contact
#' matrix will be saved as a backup in the model, automatically restoring
#' its value if it changes during the simulation.
#' @return
#' - The `set_contact_matrix()` function sets the contact matrix of the model
#' to the provided matrix and returns the modified model invisibly. The
#' function is called for its side effects and returns the modified model
#' invisibly.
set_contact_matrix <- function(model, contact_matrix, as_backup = TRUE) {
  stopifnot_model(model)
  stopifany_na(contact_matrix)
  set_contact_matrix_cpp(model, as.vector(contact_matrix), as_backup)
  invisible(model)
}
