#' Likelihood-Free Markhov Chain Monte Carlo (LFMCMC)
#'
#'
#' @aliases epiworld_lfmcmc
#' @details
#' TODO: Detail LFMCMC
#' @returns
#' - The `LFMCMC`function returns a model of class [epiworld_lfmcmc].
#' @examples
#' model_lfmcmc <- LFMCMC()
#' @export
LFMCMC <- function() {
  structure(
    LFMCMC_cpp(),
    class = c("epiworld_lfmcmc")
  )
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param params_init_ Initial model parameters
#' @param n_samples_ Number of samples
#' @param epsilon_ Epsilon parameter
#' @returns The simulated model of class `epiworld_lfmcmc`.
#' @export
run_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_) UseMethod("run_lfmcmc")

#' @export
run_lfmcmc.epiworld_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_) {
  run_lfmcmc_cpp(lfmcmc, params_init_, n_samples_, epsilon_)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param observed_data_ Observed data
#' @returns The lfmcmc model with the observed data added
#' @export
set_observed_data <- function(lfmcmc, observed_data_) UseMethod("set_observed_data")

#' @export
set_observed_data.epiworld_lfmcmc <- function(lfmcmc, observed_data_) {
  set_observed_data_cpp(lfmcmc, observed_data_)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC proposal function
#' @returns The lfmcmc model with the proposal function added
#' @export
set_proposal_fun <- function(lfmcmc, fun) UseMethod("set_proposal_fun")

#' @export
set_proposal_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_proposal_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC simulation function
#' @returns The lfmcmc model with the simulation function added
#' @export
set_simulation_fun <- function(lfmcmc, fun) UseMethod("set_simulation_fun")

#' @export
set_simulation_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_simulation_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC sumamry function
#' @returns The lfmcmc model with the summary function added
#' @export
set_summary_fun <- function(lfmcmc, fun) UseMethod("set_summary_fun")

#' @export
set_summary_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_summary_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC kernel function
#' @returns The lfmcmc model with the kernel function added
#' @export
set_kernel_fun <- function(lfmcmc, fun) UseMethod("set_kernel_fun")

#' @export
set_kernel_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_kernel_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param s The rand engine seed
#' @returns The lfmcmc model with the seed set
#' @export
seed_lfmcmc <- function(lfmcmc, s) UseMethod("seed_lfmcmc")

#' @export
seed_lfmcmc.epiworld_lfmcmc <- function(lfmcmc, s) {
  seed_lfmcmc_cpp(lfmcmc, s)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param names The model parameter names
#' @returns The lfmcmc model with the parameter names added
#' @export
set_par_names <- function(lfmcmc, names) UseMethod("set_par_names")

#' @export
set_par_names.epiworld_lfmcmc <- function(lfmcmc, names) {
  set_par_names_cpp(lfmcmc, names)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param names The model stats names
#' @returns The lfmcmc model with the stats names added
#' @export
set_stats_names <- function(lfmcmc, names) UseMethod("set_stats_names")

#' @export
set_stats_names.epiworld_lfmcmc <- function(lfmcmc, names) {
  set_stats_names_cpp(lfmcmc, names)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param x LFMCMC model to print
#' @param ... Ignored
#' @returns The lfmcmc model
#' @export
print.epiworld_lfmcmc <- function(x, ...) {
  print_lfmcmc_cpp(x)
  invisible(x)
}
