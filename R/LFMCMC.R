#' Likelihood-Free Markhov Chain Monte Carlo (LFMCMC)
#'
#'
#' @aliases epiworld_lfmcmc
#' @details
#' Performs a Likelihood-Free Markhov Chain Monte Carlo simulation
#' @param model A model of class [epiworld_model]
#' @returns
#' The `LFMCMC` function returns a model of class [epiworld_lfmcmc].
#' @examples
#' ## Setup an SIR model to use in the simulation
#' model_seed <- 122
#' model_sir <- ModelSIR(name = "COVID-19", prevalence = .1,
#'   transmission_rate = .9, recovery_rate = .3)
#' agents_smallworld(
#'   model_sir,
#'   n = 1000,
#'   k = 5,
#'   d = FALSE,
#'   p = 0.01
#' )
#' verbose_off(model_sir)
#' run(model_sir, ndays = 50, seed = model_seed)
#'
#' ## Setup LFMCMC
#' # Extract the observed data from the model
#' obs_data <- unname(as.integer(get_today_total(model_sir)))
#'
#' # Define the simulation function
#' simfun <- function(params) {
#'   set_param(model_sir, "Recovery rate", params[1])
#'   set_param(model_sir, "Transmission rate", params[2])
#'   run(model_sir, ndays = 50)
#'   res <- unname(as.integer(get_today_total(model_sir)))
#'   return(res)
#' }
#'
#' # Define the summary function
#' sumfun <- function(dat) {
#'   return(dat)
#' }
#'
#' # Create the LFMCMC model
#' lfmcmc_model <- LFMCMC(model_sir) |>
#'   set_simulation_fun(simfun) |>
#'   set_summary_fun(sumfun) |>
#'   use_proposal_norm_reflective() |>
#'   use_kernel_fun_gaussian() |>
#'   set_observed_data(obs_data)
#'
#' ## Run LFMCMC simulation
#' # Set initial parameters
#' par0 <- as.double(c(0.1, 0.5))
#' n_samp <- 2000
#' epsil <- as.double(1.0)
#'
#' # Run the LFMCMC simulation
#' run_lfmcmc(
#'   lfmcmc = lfmcmc_model,
#'   params_init_ = par0,
#'   n_samples_ = n_samp,
#'   epsilon_ = epsil,
#'   seed = model_seed
#' )
#'
#' # Print the results
#' set_stats_names(lfmcmc_model, get_states(model_sir))
#' set_par_names(lfmcmc_model, c("Immune recovery", "Infectiousness"))
#'
#' print(lfmcmc_model)
#'
#' get_stats_mean(lfmcmc_model)
#' get_params_mean(lfmcmc_model)
#'
#' @export
LFMCMC <- function(model) {
  if (!inherits(model, "epiworld_model"))
    stop("model should be of class 'epiworld_model'. It is of class ", class(model))

  structure(
    LFMCMC_cpp(model),
    class = c("epiworld_lfmcmc")
  )
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param params_init_ Initial model parameters
#' @param n_samples_ Number of samples
#' @param epsilon_ Epsilon parameter
#' @param seed Random engine seed
#' @returns The simulated model of class [epiworld_lfmcmc].
#' @export
run_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_, seed = NULL) UseMethod("run_lfmcmc")

#' @export
run_lfmcmc.epiworld_lfmcmc <- function(lfmcmc, params_init_, n_samples_, epsilon_, seed = NULL) {
  if (length(seed)) set.seed(seed)
  run_lfmcmc_cpp(lfmcmc, params_init_, n_samples_, epsilon_, sample.int(1e4, 1))
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
#' @param lfmcmc The LFMCMC model
#' @returns The LFMCMC model with proposal function set to norm reflective
#' @export
use_proposal_norm_reflective <- function(lfmcmc) {
  use_proposal_norm_reflective_cpp(lfmcmc)
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
#' @param lfmcmc The LFMCMC model
#' @returns The LFMCMC model with kernel function set to gaussian
#' @export
use_kernel_fun_gaussian <- function(lfmcmc) {
  use_kernel_fun_gaussian_cpp(lfmcmc)
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
#' @param lfmcmc LFMCMC model
#' @returns The param means for the given lfmcmc model
#' @export
get_params_mean <- function(lfmcmc) UseMethod("get_params_mean")

#' @export
get_params_mean.epiworld_lfmcmc <- function(lfmcmc) {
  get_params_mean_cpp(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @returns The stats means for the given lfmcmc model
#' @export
get_stats_mean <- function(lfmcmc) UseMethod("get_stats_mean")

#' @export
get_stats_mean.epiworld_lfmcmc <- function(lfmcmc) {
  get_stats_mean_cpp(lfmcmc)
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
