# Test just this file: tinytest::run_test_file("inst/tinytest/test-lfmcmc.R")

# Set model parameters
model_seed <- 122

# Create and run SIR Model for LFMCMC simulation -------------------------------
model_sir <- ModelSIR(name = "COVID-19", prevalence = .1, 
                      transmission_rate = .9, recovery_rate = .3)
agents_smallworld(model_sir, n = 1000, k = 5, d = FALSE, p = 0.01)
verbose_off(model_sir)
run(model_sir, ndays = 50, seed = model_seed)

# Check bad init of LFMCMC model -----------------------------------------------
expect_error(lfmcmc_bad <- LFMCMC(), 'argument "model" is missing')
expect_error(lfmcmc_bad <- LFMCMC(c("not_a_model")), "model should be of class 'epiworld_model'")

# Create LFMCMC model ----------------------------------------------------------
expect_silent(lfmcmc_model <- LFMCMC(model_sir))

# Check initialization
expect_inherits(lfmcmc_model, "epiworld_lfmcmc")
expect_length(class(lfmcmc_model), 1)

# Extract observed data from the model
obs_data <- unname(as.integer(get_today_total(model_sir)))

expect_silent(set_observed_data(lfmcmc_model, obs_data))

# Define LFMCMC functions
simfun <- function(params) {
  set_param(model_sir, "Recovery rate", params[1])
  set_param(model_sir, "Transmission rate", params[2])
  run(model_sir, ndays = 50)
  res <- unname(as.integer(get_today_total(model_sir)))
  return(res)
}

sumfun <- function(dat) { return(dat) }

propfun <- function(params_prev) {
  res <- params_prev + rnorm(length(params_prev), )
  return(res)
}

kernelfun <- function(stats_now, stats_obs, epsilon) {
  ans <- sum(mapply(function(v1, v2) (v1 - v2)^2, stats_obs, stats_now))
  return(ifelse(sqrt(ans) < epsilon, 1.0, 0.0))
}

# Check adding functions to LFMCMC
expect_silent(set_simulation_fun(lfmcmc_model, simfun))
expect_silent(set_summary_fun(lfmcmc_model, sumfun))
expect_silent(set_proposal_fun(lfmcmc_model, propfun))
expect_silent(set_kernel_fun(lfmcmc_model, kernelfun))

# Run LFMCMC simulation --------------------------------------------------------
# Initial parameters
par0 <- as.double(c(0.1, 0.5))
n_samp <- 2000
epsil <- as.double(1.0)

expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  n_samples_ = n_samp,
  epsilon_ = epsil,
  seed = model_seed
))

expect_silent(set_stats_names(lfmcmc_model, get_states(model_sir)))
expect_silent(set_par_names(lfmcmc_model, c("Immune recovery", "Infectiousness")))

expect_stdout(print(lfmcmc_model))

expect_equal(get_stats_mean(lfmcmc_model), c(4.45, 2.6135, 992.4365))
expect_equal(get_params_mean(lfmcmc_model), c(11.58421, 18.96851), tolerance = 0.00001)

# Check LFMCMC using factory functions -----------------------------------------
expect_silent(use_proposal_norm_reflective(lfmcmc_model))
expect_silent(use_kernel_fun_gaussian(lfmcmc_model))

expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  n_samples_ = n_samp,
  epsilon_ = epsil,
  seed = model_seed
))

# Check running LFMCMC with missing parameters ---------------------------------
expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  n_samples_ = n_samp,
  epsilon_ = epsil
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  n_samples_ = n_samp
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  epsilon_ = epsil
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  n_samples_ = n_samp,
  epsilon_ = epsil
))

expect_error(run_lfmcmc(
  params_init_ = par0,
  n_samples_ = n_samp,
  epsilon_ = epsil
))

expect_error(run_lfmcmc(lfmcmc = lfmcmc_model))
