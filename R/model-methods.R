#' Common methods for predefined models of Epiworld
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
#' @rdname epiworld-methods
agents_smallworld <- function(m, n, k, d, p) UseMethod("agents_smallworld")

#' @export
#' @rdname epiworld-methods
run <- function(m) UseMethod("run")
