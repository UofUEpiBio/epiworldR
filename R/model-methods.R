#' Common methods for predefined models of Epiworld
#' @export
#' @name epiworld-methods
init <- function(m, days, seed) UseMethod("init")

#' @export
#' @rdname epiworld-methods
agents_smallworld <- function(m, n, k, d, p) UseMethod("agents_smallworld")

#' @export
#' @rdname epiworld-methods
run <- function(m) UseMethod("run")
