#' Likelihood-Free Markhov Chain Monte Carlo (LFMCMC)
#'
#' @aliases epiworld_lfmcmc
#' @param model A model of class [epiworld_model] or `NULL` (see details).
#' @param fun A function (see details).
#' @details
#' Performs a Likelihood-Free Markhov Chain Monte Carlo simulation. When
#' `model` is not `NULL`, the model uses the same random-number generator
#' engine as the model. Otherwise, when `model` is `NULL`, a new random-number
#' generator engine is created.
#'
#' The functions passed to the LFMCMC object have different arguments depending
#' on the object:
#' - `set_proposal_fun`: A vector of parameters and the model.
#' - `set_simulation_fun`: A vector of parameters and the model.
#' - `set_summary_fun`: A vector of simulated data and the model.
#' - `set_kernel_fun`: A vector of simulated statistics, observed statistics,
#' epsilon, and the model.
#'
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
#' obs_data <- get_today_total(model_sir)
#'
#' # Define the simulation function
#' simfun <- function(params, lfmcmc_obj) {
#'   set_param(model_sir, "Recovery rate", params[1])
#'   set_param(model_sir, "Transmission rate", params[2])
#'   run(model_sir, ndays = 50)
#'   res <- get_today_total(model_sir)
#'   return(res)
#' }
#'
#' # Define the summary function
#' sumfun <- function(dat, lfmcmc_obj) {
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
#' par0 <- c(0.1, 0.5)
#' n_samp <- 2000
#' epsil <- 1.0
#'
#' # Run the LFMCMC simulation
#' verbose_off(lfmcmc_model)
#' run_lfmcmc(
#'   lfmcmc = lfmcmc_model,
#'   params_init = par0,
#'   n_samples = n_samp,
#'   epsilon = epsil,
#'   seed = model_seed
#' )
#'
#' # Print the results
#' set_stats_names(lfmcmc_model, get_states(model_sir))
#' set_params_names(lfmcmc_model, c("Immune recovery", "Infectiousness"))
#'
#' print(lfmcmc_model)
#'
#' get_mean_stats(lfmcmc_model)
#' get_mean_params(lfmcmc_model)
#'
#' @export
#' @concept fmcmc
LFMCMC <- function(model = NULL) {

  if ((length(model) > 0) && !inherits(model, "epiworld_model"))
    stop(
      "model should be of class 'epiworld_model'. It is of class ",
      paste(class(model), collapse = "\", ")
    )

  structure(
    LFMCMC_cpp(model),
    class = c("epiworld_lfmcmc")
  )

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param params_init Initial model parameters, treated as double
#' @param n_samples Number of samples, treated as integer
#' @param epsilon Epsilon parameter, treated as double
#' @param seed Random engine seed
#' @returns The simulated model of class [epiworld_lfmcmc].
#' @export
run_lfmcmc <- function(
    lfmcmc, params_init, n_samples, epsilon,
    seed = NULL
    ) {

  stopifnot_lfmcmc(lfmcmc)

  if (length(seed))
    set.seed(seed)

  run_lfmcmc_cpp(
    lfmcmc,
    as.double(params_init),
    as.integer(n_samples),
    as.double(epsilon),
    sample.int(1e4, 1)
  )

  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param observed_data Observed data, treated as double.
#' @export
set_observed_data <- function(lfmcmc, observed_data) {

  stopifnot_lfmcmc(lfmcmc)

  set_observed_data_cpp(
    lfmcmc,
    as.double(observed_data)
  )

  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @export
set_proposal_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_proposal_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @export
use_proposal_norm_reflective <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  use_proposal_norm_reflective_cpp(lfmcmc)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @export
set_simulation_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_simulation_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @export
set_summary_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_summary_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @export
set_kernel_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_kernel_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @returns
#' - `use_kernel_fun_gaussian`: The LFMCMC model with kernel function set to
#' gaussian.
#' @export
use_kernel_fun_gaussian <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  use_kernel_fun_gaussian_cpp(lfmcmc)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @returns
#' - `get_mean_params`: The param means for the given lfmcmc model.
#' @export
get_mean_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_mean_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @returns
#' - `get_mean_stats`: The stats means for the given lfmcmc model.
#' @export
get_mean_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_mean_stats_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_initial_params` returns the initial parameters
#' for the given LFMCMC model.
get_initial_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_initial_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_current_proposed_params` returns the proposed parameters
#' for the next LFMCMC sample.
get_current_proposed_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_current_proposed_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_current_accepted_params` returns the most recently accepted
#' parameters (the current state of the LFMCMC)
get_current_accepted_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_current_accepted_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_current_proposed_stats` returns the statistics
#' from the simulation run with the proposed parameters
get_current_proposed_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_current_proposed_stats_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_current_accepted_stats` returns the statistics
#' from the most recently accepted parameters
get_current_accepted_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_current_accepted_stats_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_observed_stats` returns the statistics
#' for the observed data
get_observed_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_observed_stats_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_sample_params` returns a matrix of sample
#' parameters for the given LFMCMC model. with the number of rows equal to the
#' number of samples and the number of columns equal to the number of
#' parameters.
get_all_sample_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  a_params <- get_all_sample_params_cpp(lfmcmc)
  n_params <- get_n_params(lfmcmc)

  matrix(
    a_params,
    ncol = n_params,
    byrow = TRUE
  )

}

#' @export
#' @rdname LFMCMC
#' @returns
#' - The function `get_all_sample_stats` returns a matrix of statistics
#' for the given LFMCMC model. with the number of rows equal to the number of
#' samples and the number of columns equal to the number of statistics.
get_all_sample_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  stats <- get_all_sample_stats_cpp(lfmcmc)
  n_stats <- get_n_stats(lfmcmc)

  matrix(
    stats,
    ncol = n_stats,
    byrow = TRUE
  )

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_sample_acceptance` returns a vector of boolean flags
#' which indicate whether a given sample was accepted
get_all_sample_acceptance <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_all_sample_acceptance_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_sample_drawn_prob` returns a vector of drawn probabilities
#' for each sample
get_all_sample_drawn_prob <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_all_sample_drawn_prob_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_sample_kernel_scores` returns a vector of kernel scores for
#' each sample
get_all_sample_kernel_scores <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_all_sample_kernel_scores_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_accepted_params` returns a matrix of accepted
#' parameters for the given LFMCMC model. with the number of rows equal to the
#' number of samples and the number of columns equal to the number of
#' parameters.
get_all_accepted_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  a_params <- get_all_accepted_params_cpp(lfmcmc)
  n_params <- get_n_params(lfmcmc)

  matrix(
    a_params,
    ncol = n_params,
    byrow = TRUE
  )

}


#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_accepted_stats` returns a matrix of accepted statistics
#' for the given LFMCMC model. with the number of rows equal to the number of
#' samples and the number of columns equal to the number of statistics.
get_all_accepted_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  a_stats <- get_all_accepted_stats_cpp(lfmcmc)
  n_stats <- get_n_stats(lfmcmc)

  matrix(
    a_stats,
    ncol = n_stats,
    byrow = TRUE
  )

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_all_accepted_kernel_scores` returns a vector of kernel scores for
#' each accepted sample
get_all_accepted_kernel_scores <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_all_accepted_kernel_scores_cpp(lfmcmc)

}

#' @export
#' @rdname LFMCMC
#' @returns
#' - The functions `get_n_samples`, `get_n_stats`, and `get_n_params`
#' return the number of samples, statistics, and parameters for the given
#' LFMCMC model, respectively.
get_n_samples <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_samples_cpp(lfmcmc)

}

#' @export
#' @rdname LFMCMC
get_n_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_stats_cpp(lfmcmc)

}

#' @export
#' @rdname LFMCMC
get_n_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The `verbose_on` and `verbose_off` functions return the same model, however
#' `verbose_off` returns the model with no progress bar.
#' @details
#' The `verbose_on` and `verbose_off` functions activate and deactivate printing
#' progress on screen, respectively. Both functions return the model (`x`) invisibly.

#' @export
verbose_off.epiworld_lfmcmc <- function(x) {
  invisible(verbose_off_lfmcmc_cpp(x))
}

#' @export
verbose_on.epiworld_lfmcmc <- function(x) {
  invisible(verbose_on_lfmcmc_cpp(x))
}

#' @rdname LFMCMC
#' @param names Character vector of names.
#' @returns
#' - `set_params_names`: The lfmcmc model with the parameter names added.
#' @export
set_params_names <- function(lfmcmc, names) {

  stopifnot_lfmcmc(lfmcmc)
  set_params_names_cpp(lfmcmc, names)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @returns
#' - `set_stats_names`: The lfmcmc model with the stats names added.
#' @export
set_stats_names <- function(lfmcmc, names) {

  stopifnot_lfmcmc(lfmcmc)
  set_stats_names_cpp(lfmcmc, names)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param x LFMCMC model to print
#' @param ... Ignored
#' @param burnin Integer. Number of samples to discard as burnin before
#' computing the summary.
#' @export
print.epiworld_lfmcmc <- function(x, burnin = 0, ...) {

  if (!is.numeric(burnin))
    stop("The 'burnin' argument must be an integer.")

  if (burnin < 0)
    stop("The 'burnin' argument must be a non-negative integer.")

  print_lfmcmc_cpp(x, burnin = burnin)
  invisible(x)

}
