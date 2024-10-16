#' LFMCMC
#'
#'
#' @export
LFMCMC <- function() {
  structure(
    LFMCMC_cpp(),
    class = c("epiworld_lfmcmc")
  )
}

#' @export
run.epiworld_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_) {
  run_lfmcmc_cpp(lfmcmc, params_init_, n_samples_, epsilon_)
  invisible(lfmcmc)
}

# #' @export
# run.epiworld_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_) {
#   run_lfmcmc_cpp(lfmcmc, params_init_, n_samples_, epsilon_)
#   invisible(lfmcmc)
# }

# #' @export
# set_observed_data.epiworld_lfmcmc <- function(lfmcmc, observed_data_) {
#   set_observed_data_cpp(lfmcmc, observed_data_)
#   invisible(lfmcmc)
# }

# #' @export
# set_proposal_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
#   set_proposal_fun_cpp(lfmcmc, fun)
#   invisible(lfmcmc)
# }

# #' @export
# set_simulation_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
#   set_simulation_fun_cpp(lfmcmc, fun)
#   invisible(lfmcmc)
# }

# #' @export
# set_summary_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
#   set_summary_fun_cpp(lfmcmc, fun)
#   invisible(lfmcmc)
# }

# #' @export
# set_kernel_fun.epiworld_lfmcmc <- function(lfmcmc, fun) {
#   set_kernel_fun_cpp(lfmcmc, fun)
#   invisible(lfmcmc)
# }
