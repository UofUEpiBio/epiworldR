# Test just this file: tinytest::run_test_file("inst/tinytest/test-lfmcmc.R")

# Set model parameters
model_seed <- 122

# Create and run SIR Model for LFMCMC simulation -------------------------------
model_sir <- ModelSIR(name = "COVID-19", prevalence = .1, 
                      transmission_rate = .9, recovery_rate = .3)
agents_smallworld(model_sir, n = 1000, k = 5, d = FALSE, p = 0.01)
verbose_off(model_sir)
run(model_sir, ndays = 50, seed = model_seed)

# Create LFMCMC model ----------------------------------------------------------
lfmcmc_model <- LFMCMC(model_sir)

# Check initialization
expect_inherits(lfmcmc_model, "epiworld_lfmcmc")

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

# Create LFMCMC simulation -----------------------------------------------------
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