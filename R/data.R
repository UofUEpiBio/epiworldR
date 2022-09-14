#' Accessing the database of `epiworld`
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' any model.
#' @name epiworld-data
#' @family Models
NULL

#' @export
#' @rdname epiworld-data
get_hist_total <- function(x) get_hist_total_cpp(x, class(x))

#' @export
#' @rdname epiworld-data
get_transition_probability <- function(x) {
  res <- get_transition_probability_cpp(x, class(x))
  s   <- get_status(x)
  
  ns <- length(s)
  
  matrix(res, nrow = ns, ncol = ns, dimnames = list(s, s))
}

#' @export
#' @rdname epiworld-data
get_status <- function(x) get_status_cpp(x, class(x))

#' @export
#' @rdname epiworld-data
get_reproductive_number <- function(x) get_reproductive_number_cpp(x, class(x))

