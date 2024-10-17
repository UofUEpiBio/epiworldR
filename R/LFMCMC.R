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
run_lfmcmc.epiworld_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_) {
  run_lfmcmc_cpp(lfmcmc, params_init_, n_samples_, epsilon_)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param observed_data_ Observed data
#' @returns The lfmcmc model with the observed data added
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
set_proposal_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_proposal_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC simulation function
#' @returns The lfmcmc model with the simulation function added
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
set_summary_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_summary_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC kernel function
#' @returns The lfmcmc model with the kernel function added
#' @export
set_kernel_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
  set_kernel_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)
}
