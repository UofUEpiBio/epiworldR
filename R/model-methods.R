#' Common methods for predefined models of Epiworld
#' @param x An object of class `epiworld_model`.
#' @param ndays Number of days (steps) of the simulation.
#' @param seed Seed to set for initializing random number generator.
#' @param m Model object.
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

#' @name epiworld-methods
#' @export
#' @details
#' The `verbose_on` and `verbose_off` functions activate and deactivate printing
#' progress on screen, respectively. Both functions return the model (`x`) invisibly.
verbose_off <- function(x) UseMethod("verbose_off")

#' @export
verbose_off.epiworld_model <- function(x) {
  verbose_off_cpp(x)
  invisible(x)
}

#' @name epiworld-methods
#' @export
verbose_on <- function(x) UseMethod("verbose_on")

#' @export
verbose_on.epiworld_model <- function(x) {
  verbose_on_cpp(x)
  invisible(x)
}

#' @export
#' @rdname epiworld-methods
run <- function(m, ndays, seed = sample.int(1e4, 1)) UseMethod("run")

#' @export
run.epiworld_model <- function(m, ndays, seed = sample.int(1e4, 1)) {
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
get_states <- function(x) UseMethod("get_states")

#' @export
get_states.epiworld_model <- function(x) get_states_cpp(x)

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
#' @returns
#' - `get_name` returns the name of the model.
#' @rdname epiworld-methods
get_name <- function(x) UseMethod("get_name")

#' @export
get_name.epiworld_model <- function(x) {
  get_name_cpp(x)
}

#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_n_variants` returns the number of variants of the model.
get_n_variants <- function(x) UseMethod("get_n_variants")

#' @export
get_n_variants.epiworld_model <- function(x) get_n_variants_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_n_tools` returns the number of tools of the model.
get_n_tools <- function(x) UseMethod("get_n_tools")

#' @export
get_n_tools.epiworld_model <- function(x) get_n_tools_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_ndays` returns the number of days of the model.
get_ndays <- function(x) UseMethod("get_ndays")

#' @export
get_ndays.epiworld_model <- function(x) get_ndays_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns 
#' - `get_n_replicates` returns the number of replicates of the model.
get_n_replicates <- function(x) UseMethod("get_n_replicates")

#' @export
get_n_replicates.epiworld_model <- function(x) get_n_replicates_cpp(x)


#' @export
#' @rdname epiworld-methods
#' @returns 
#' - `size.epiworld_model` returns the number of agents in the model.
#' 
size <- function(x) UseMethod("size")

#' @export
size.epiworld_model <- function(x) size_cpp(x)

