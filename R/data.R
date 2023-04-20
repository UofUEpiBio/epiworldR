#' Accessing the database of `epiworld`
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' any model.
#' @param skip_zeros Logical scalar. When `FALSE` it will return all the
#' entries in the transition matrix.
#' @param ... In the case of `plot.epiworld_repnum`, further arguments passed to
#' [graphics::plot].
#' @name epiworld-data
#' @family Models
#' @examples
#' # SEIR Connected
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0, 
#'   prob_transmission = 0.8,
#'   incubation_days   = 7.0,
#'   prob_recovery     = 0.3
#' )
#' 
#' # Running the simulation for 50 steps (days)
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Retrieving the transition probability
#' get_transition_probability(seirconn) 
#' 
#' # We can go further and get all the history
#' t_hist <- get_hist_transition_matrix(seirconn)
#' 
#' head(t_hist)
#' 
#' # And turn it into an array
#' as.array(t_hist)[,,1:3]
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
get_hist_variant <- function(x) UseMethod("get_hist_variant")

#' @export
get_hist_variant.epiworld_model <- function(x)  {
  get_hist_variant_cpp(x)
}

#' @export
#' @rdname epiworld-data
get_hist_tool <- function(x) UseMethod("get_hist_tool")

#' @export
get_hist_tool.epiworld_model <- function(x)  {
  get_hist_tool_cpp(x)
}

#' @export
#' @rdname epiworld-data
get_transition_probability <- function(x) {
  UseMethod("get_transition_probability")
}

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
get_reproductive_number.epiworld_model <- function(x) {
  res <- get_reproductive_number_cpp(x)
  class(res) <- c("epiworld_repnum", class(res))
  res
}

#' @rdname epiworld-data
#' @export
plot.epiworld_repnum <- function(
    x,
    y = NULL,
    ylab = "Average Rep. Number",
    xlab = "Day (step)",
    main = "Daily Average Reproductive Number",
    ...) {
  
  res <- tapply(
    X     = x[, "counts"],
    INDEX = x[, "source_exposure_dates"],
    FUN   = mean
    )
  
  graphics::plot(
    x = sort(as.integer(names(res))),
    y = res,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
    )
  
}

#' @export
#' @rdname epiworld-data
#' @returns 
#' In the case of `get_hist_transition_matrix`, a [data.frame] with four columns:
#' "state_from", "state_to", "date", and "counts." It will be
get_hist_transition_matrix <- function(x, skip_zeros = FALSE)
  UseMethod("get_hist_transition_matrix")

#' @export
get_hist_transition_matrix.epiworld_model <- function(x, skip_zeros = FALSE) {
  
  res <- get_hist_transition_matrix_cpp(x, skip_zeros)
  class(res) <- c(class(res), "epiworld_hist_transition")
  
  attr(res, "states") <- get_state(x)
  attr(res, "nsteps") <- get_ndays(x)
  
  res
  
}


#' @export
#' @returns 
#' The `as.array` method for `epiworld_hist_transition` objects turns the `data.frame`
#' returned by `get_hist_transition_matrix` into an array of `nstates x nstates x (ndays + 1)`
#' entries, where the first entry is the initial state.
#' @rdname epiworld-data
as.array.epiworld_hist_transition <- function(x, ...) {
  
  states <- attr(x, "states")
  n_states <- length(states)
  n_steps  <- attr(x, "nsteps")
  
  res <- array(
    0L,
    dim      = c(n_states, n_states, n_steps + 1), # Includes the baseline
    dimnames = list(states, states, 0:n_steps)
    )
 
  res[cbind(x[,1], x[,2], x[,3])] <- x[,4]
  
  res

}

