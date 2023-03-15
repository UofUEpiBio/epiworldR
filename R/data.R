#' Accessing the database of `epiworld`
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' any model.
#' @name epiworld-data
#' @family Models
NULL

#' @export
#' @rdname epiworld-data
get_hist_total <- function(x) UseMethod("get_hist_total")

#' @export
get_hist_total.epiworld_model <- function(x)  {
  get_hist_total_cpp(x)
}

#' @export
#' @rdname epiworld-data
get_transition_probability <- function(x) UseMethod("get_transition_probability")

#' @export
get_transition_probability.epiworld_model <- function(x)  {
  res <- get_transition_probability_cpp(x)
  s   <- get_state(x)
  
  ns <- length(s)
  
  matrix(res, nrow = ns, ncol = ns, dimnames = list(s, s))
}

#' @export
#' @rdname epiworld-data
get_reproductive_number <- function(x) UseMethod("get_reproductive_number")

#' @export
get_reproductive_number.epiworld_model <- function(x) get_reproductive_number_cpp(x)

