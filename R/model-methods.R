#' Common methods for predefined models of Epiworld
#' @param x An object of class `epiworld_model`.
#' @param days Number of days.
#' @param seed Seed to set for initializing random number generator.
#' @param m Model object.
#' @param n Number of individuals in the population.
#' @param k Number of ties in the small world network.
#' @param d Logical scalar. Whether the graph is directed or not.
#' @param p Probability of rewiring.
#' @export
#' @name epiworld-methods
init <- function(m, days, seed) UseMethod("init")

#' @export
init.epiworld_model <- function(m, days, seed) {
  init_cpp(m, days, seed)
  invisible(m)
}


#' @name epiworld-methods
#' @export
queuing_on <- function(x) UseMethod("queuing_on")

#' @export
queuing_on.epiworld_model <- function(x) {
  queuing_on_cpp(x)
  invisible(x)
}

#' @name epiworld-methods
#' @export
queuing_off <- function(x) UseMethod("queuing_off")

#' @export
queuing_off.epiworld_model <- function(x) {
  queuing_off_cpp(x)
  invisible(x)
}

#' @export
#' @rdname epiworld-methods
agents_smallworld <- function(m, n, k, d, p) UseMethod("agents_smallworld")

#' @export
agents_smallworld.epiworld_model <- function(m, n, k, d, p) {
  agents_smallworld_cpp(m, n, k, d, p)
  invisible(m)
}

#' @export
#' @rdname epiworld-methods
run <- function(m) UseMethod("run")

#' @export
run.epiworld_model <- function(m) {
  run_cpp(m)
  invisible(m)
}

#' @export
print.epiworld_model <- function(x, ...) {
  print_cpp(x)
  invisible(x)
}



