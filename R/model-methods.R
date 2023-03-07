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
run <- function(m, ndays, seed = sample.int(1e4, 1)) UseMethod("run")

#' @export
run.epiworld_model <- function(m, ndays, seed) {
  run_cpp(m, ndays, seed)
  invisible(m)
}

#' @export
print.epiworld_model <- function(x, ...) {
  print_cpp(x)
  invisible(x)
}

#' @export
#' @rdname epiworld-methods
get_status <- function(x) UseMethod("get_status")

#' @export
get_status.epiworld_model <- function(x) get_status_cpp(x)

#' @export
#' @param pname String, name of the parameter
#' @rdname epiworld-methods
get_param <- function(x, pname) UseMethod("get_param")

#' @export
get_param.epiworld_model <- function(x, pname) {
  get_param_cpp(x, pname)
}


#' @export
#' @param pval Numeric. Value of the parameter
#' @rdname epiworld-methods
set_param <- function(x, pname, pval) UseMethod("set_param")

#' @export
set_param.epiworld_model <- function(x, pname, pval) {
  invisible(set_param_cpp(x, pname, pval))
  invisible(x)
}

#' @export
#' @param mname String. Name of the model
#' @rdname epiworld-methods
set_name <- function(x, mname) UseMethod("set_name")

#' @export
set_name.epiworld_model <- function(x, mname) {
  set_name_cpp(x, mname)
  invisible(x)
}

#' @export
#' @rdname epiworld-methods
get_name <- function(x) UseMethod("get_name")

#' @export
get_name.epiworld_model <- function(x) {
  get_name_cpp(x)
}

