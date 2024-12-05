stopifnot_lfmcmc <- function(x) {
  # Catching the value of x
  nam <- match.call()$x

  if (!inherits(x, "epiworld_lfmcmc"))
    stop(nam, " must be an object of class epiworld_lfmcmc")

}


#' Likelihood-Free Markhov Chain Monte Carlo (LFMCMC)
#'
#' @aliases epiworld_lfmcmc
#' @param model A model of class [epiworld_model] or `NULL` (see details).
#' @details
#' Performs a Likelihood-Free Markhov Chain Monte Carlo simulation. When
#' `model` is not `NULL`, the model uses the same random-number generator
#' engine as the model. Otherwise, when `model` is `NULL`, a new random-number
#' generator engine is created.
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
#' simfun <- function(params) {
#'   set_param(model_sir, "Recovery rate", params[1])
#'   set_param(model_sir, "Transmission rate", params[2])
#'   run(model_sir, ndays = 50)
#'   res <- get_today_total(model_sir)
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
#' par0 <- c(0.1, 0.5)
#' n_samp <- 2000
#' epsil <- 1.0
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
#' set_stat_names(lfmcmc_model, get_states(model_sir))
#' set_param_names(lfmcmc_model, c("Immune recovery", "Infectiousness"))
#'
#' print(lfmcmc_model)
#'
#' get_mean_stats(lfmcmc_model)
#' get_mean_params(lfmcmc_model)
#'
#' @export
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
#' @param params_init_ Initial model parameters, treated as double
#' @param n_samples_ Number of samples, treated as integer
#' @param epsilon_ Epsilon parameter, treated as double
#' @param seed Random engine seed
#' @returns The simulated model of class [epiworld_lfmcmc].
#' @export
run_lfmcmc <- function(
    lfmcmc, params_init_, n_samples_, epsilon_,
    seed = NULL
    ) {

  stopifnot_lfmcmc(lfmcmc)

  if (length(seed))
    set.seed(seed)

  run_lfmcmc_cpp(
    lfmcmc,
    as.double(params_init_),
    as.integer(n_samples_),
    as.double(epsilon_),
    sample.int(1e4, 1)
  )

  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param observed_data_ Observed data, treated as double
#' @returns The lfmcmc model with the observed data added
#' @export
set_observed_data <- function(lfmcmc, observed_data_) {

  stopifnot_lfmcmc(lfmcmc)

  set_observed_data_cpp(
    lfmcmc,
    as.double(observed_data_)
  )

  invisible(lfmcmc)
}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC proposal function
#' @export
#' @returns The lfmcmc model with the proposal function added
set_proposal_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_proposal_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc The LFMCMC model
#' @returns The LFMCMC model with proposal function set to norm reflective
#' @export
use_proposal_norm_reflective <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  use_proposal_norm_reflective_cpp(lfmcmc)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC simulation function
#' @returns The lfmcmc model with the simulation function added
#' @export
set_simulation_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_simulation_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC sumamry function
#' @returns The lfmcmc model with the summary function added
#' @export
set_summary_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_summary_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param fun The LFMCMC kernel function
#' @returns The lfmcmc model with the kernel function added
#' @export
set_kernel_fun <- function(lfmcmc, fun) {

  stopifnot_lfmcmc(lfmcmc)
  set_kernel_fun_cpp(lfmcmc, fun)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc The LFMCMC model
#' @returns The LFMCMC model with kernel function set to gaussian
#' @export
use_kernel_fun_gaussian <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  use_kernel_fun_gaussian_cpp(lfmcmc)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param names The model parameter names
#' @returns The lfmcmc model with the parameter names added
#' @export
set_param_names <- function(lfmcmc, names) {

  stopifnot_lfmcmc(lfmcmc)
  set_param_names_cpp(lfmcmc, names)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @param names The model stats names
#' @returns The lfmcmc model with the stats names added
#' @export
set_stat_names <- function(lfmcmc, names) {

  stopifnot_lfmcmc(lfmcmc)
  set_stat_names_cpp(lfmcmc, names)
  invisible(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @returns The param means for the given lfmcmc model
#' @export
get_mean_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_mean_params_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @param lfmcmc LFMCMC model
#' @returns The stats means for the given lfmcmc model
#' @export
get_mean_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_mean_stats_cpp(lfmcmc)

}

#' @rdname LFMCMC
#' @param x LFMCMC model to print
#' @param ... Ignored
#' @param burnin Integer. Number of samples to discard as burnin before computing the summary.
#' @returns The lfmcmc model
#' @export
print.epiworld_lfmcmc <- function(x, burnin = 0, ...) {

  if (!is.numeric(burnin))
    stop("The 'burnin' argument must be an integer.")

  if (burnin < 0)
    stop("The 'burnin' argument must be a non-negative integer.")

  print_lfmcmc_cpp(x, burnin = burnin)
  invisible(x)

}

#' @rdname LFMCMC
#' @export
#' @returns
#' - The function `get_accepted_params` returns a matrix of accepted
#' parameters for the given LFMCMC model. with the number of rows equal to the
#' number of samples and the number of columns equal to the number of
#' parameters.
get_accepted_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  a_params <- get_accepted_params_cpp(lfmcmc)
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
#' - The function `get_accepted_stats` returns a matrix of accepted statistics
#' for the given LFMCMC model. with the number of rows equal to the number of
#' samples and the number of columns equal to the number of statistics.
get_accepted_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  a_stats <- get_accepted_stats_cpp(lfmcmc)
  n_stats <- get_n_stats(lfmcmc)

  matrix(
    a_stats,
    ncol = n_stats,
    byrow = TRUE
  )

}

#' @export
#' @rdname LFMCMC
#' @returns
#' - The function `get_sample_stats` returns a matrix of statistics
#' for the given LFMCMC model. with the number of rows equal to the number of
#' samples and the number of columns equal to the number of statistics.
get_sample_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  stats <- get_sample_stats_cpp(lfmcmc)
  n_stats <- get_n_stats(lfmcmc)

  matrix(
    stats,
    ncol = n_stats,
    byrow = TRUE
  )

}

#' @export
#' @rdname LFMCMC
#' @returns
#' - The functions `get_n_params`, `get_n_stats`, and `get_n_samples`
#' return the number of parameters, statistics, and samples for the given
#' LFMCMC model, respectively.
get_n_params <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_params_cpp(lfmcmc)

}

#' @export
#' @rdname LFMCMC
get_n_stats <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_stats_cpp(lfmcmc)

}

#' @export
#' @rdname LFMCMC
get_n_samples <- function(lfmcmc) {

  stopifnot_lfmcmc(lfmcmc)
  get_n_samples_cpp(lfmcmc)

}
