# Test just this file: tinytest::run_test_file("inst/tinytest/test-lfmcmc.R")

# Set model parameters
model_seed <- 122

# Create and run SIR Model for LFMCMC simulation -------------------------------
model_sir <- ModelSIR(name = "COVID-19", prevalence = .1, 
                      transmission_rate = .3, recovery_rate = .3)
agents_smallworld(model_sir, n = 1000, k = 5, d = FALSE, p = 0.01)
verbose_off(model_sir)
run(model_sir, ndays = 50, seed = model_seed)

# Check init of LFMCMC model without epiworld model ----------------------------
expect_silent(lfmcmc_nomodel <- LFMCMC())

# Check bad init of LFMCMC model -----------------------------------------------
expect_error(lfmcmc_bad <- LFMCMC(c("not_a_model")), "model should be of class 'epiworld_model'")

# Create LFMCMC model ----------------------------------------------------------
expect_silent(lfmcmc_model <- LFMCMC(model_sir))

# Check initialization
expect_inherits(lfmcmc_model, "epiworld_lfmcmc")
expect_length(class(lfmcmc_model), 1)

# Extract observed data from the model
obs_data <- get_today_total(model_sir)

expect_silent(set_observed_data(lfmcmc_model, obs_data))

# Define LFMCMC functions
simfun <- function(params, model) {
  set_param(model_sir, "Recovery rate", params[1])
  set_param(model_sir, "Transmission rate", params[2])
  run(model_sir, ndays = 50)
  res <- get_today_total(model_sir)
  return(res)
}

sumfun <- function(dat, model) { return(dat) }

propfun <- function(old_params, model) {
  res <- plogis(qlogis(old_params) + rnorm(length(old_params)))
  return(res)
}

kernelfun <- function(simulated_stats, observed_stats, epsilon, model) {
  dnorm(sqrt(sum((simulated_stats - observed_stats)^2)))
}

# Check adding functions to LFMCMC
expect_silent(set_simulation_fun(lfmcmc_model, simfun))
expect_silent(set_summary_fun(lfmcmc_model, sumfun))
expect_silent(set_proposal_fun(lfmcmc_model, propfun))
expect_silent(set_kernel_fun(lfmcmc_model, kernelfun))

# Run LFMCMC simulation --------------------------------------------------------
# Initial parameters
par0 <- c(0.1, 0.5)
n_samp <- 2000
epsil <- 1.0

expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0,
  n_samples = n_samp,
  epsilon = epsil,
  seed = model_seed
))

expect_silent(set_stats_names(lfmcmc_model, get_states(model_sir)))
expect_silent(set_params_names(lfmcmc_model, c("Immune recovery", "Infectiousness")))

# Check printing LFMCMC --------------------------------------------------------
expect_stdout(print(lfmcmc_model))
expect_stdout(print(lfmcmc_model, burnin = n_samp / 2))
expect_error(print(lfmcmc_model, burnin = n_samp), "burnin is greater than or equal to the number of samples")
expect_error(print(lfmcmc_model, burnin = n_samp + 50), "burnin is greater than or equal to the number of samples")
expect_error(print(lfmcmc_model, burnin = -n_samp / 2), "argument must be a non-negative integer")
expect_error(print(lfmcmc_model, burnin = "n_samp"), "argument must be an integer")

# Check LFMCMC getters ---------------------------------------------------------
expect_equal(get_n_samples(lfmcmc_model), n_samp)

expected_stats_mean <- c(284.7140, 0.8485, 713.9375)
expect_equal(get_mean_stats(lfmcmc_model), expected_stats_mean)
expect_equal(get_n_stats(lfmcmc_model), length(expected_stats_mean))

expected_params_mean <- c(0.3133401, 0.2749686)
expect_equal(get_mean_params(lfmcmc_model), expected_params_mean, tolerance = 0.0001)
expect_equal(get_n_params(lfmcmc_model), length(expected_params_mean))

expect_equal(dim(get_accepted_params(lfmcmc_model)), c(n_samp, length(expected_params_mean)))
expect_equal(dim(get_sample_stats(lfmcmc_model)), c(n_samp, length(expected_stats_mean)))

# Check LFMCMC using factory functions -----------------------------------------
expect_silent(use_proposal_norm_reflective(lfmcmc_model))
expect_silent(use_kernel_fun_gaussian(lfmcmc_model))

expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0,
  n_samples = n_samp,
  epsilon = epsil,
  seed = model_seed
))

# Check LFMCMC type coercion of parameters and observed data -------------------
obs_data_int <- as.integer(obs_data)
expect_silent(set_observed_data(lfmcmc_model, obs_data_int))

par0_int <- as.integer(c(1, 5))
n_samp_double <- as.double(2000.0)
epsil_int <- as.integer(1)

expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0_int,
  n_samples = n_samp_double,
  epsilon = epsil_int,
  seed = model_seed
))

# Check running LFMCMC with missing parameters ---------------------------------
expect_silent(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0,
  n_samples = n_samp,
  epsilon = epsil
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0,
  n_samples = n_samp
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = par0,
  epsilon = epsil
))

expect_error(run_lfmcmc(
  lfmcmc = lfmcmc_model,
  n_samples = n_samp,
  epsilon = epsil
))

expect_error(run_lfmcmc(
  params_init = par0,
  n_samples = n_samp,
  epsilon = epsil
))

expect_error(run_lfmcmc(lfmcmc = lfmcmc_model))

# Check running LFMCMC without epiworld model ----------------------------------
model_seed <- 4222
set.seed(model_seed)
Y <- rnorm(2000, mean = -5, sd = 2.5)

# Define LFMCMC functions
simfun <- function(par, model) {
  rnorm(2000, mean = par[1], sd = par[2])
}

sumfun <- function(x, model) {
  c(mean(x), sd(x))
}

propfun <- function(par, model) {
  
  par_new <- par + rnorm(2, sd = 0.1)

  # Reflecting par2
  if (par_new[2] < 0) {
    par_new[2] <- par[2] - (par_new[2] - par[2])
  }

  return(par_new)
}

kernelfun <- function(simulated_stats, observed_stats, epsilon, model) {

  dnorm(sqrt(sum((observed_stats - simulated_stats)^2)))
  
}

# Setup LFMCMC
lfmcmc_model <- LFMCMC() |>
  set_simulation_fun(simfun) |>
  set_summary_fun(sumfun) |>
  set_proposal_fun(propfun) |>
  set_kernel_fun(kernelfun) |>
  set_observed_data(Y)

# Run LFMCMC
x <- run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init = c(0, 1),
  n_samples = 3000,
  epsilon = 1.0,
  seed = model_seed
)

x_means <- get_mean_params(x)

expect_equivalent(
  x_means,
  c(-5, 2.5),
  tolerance = 0.1
)

# Check functions fail when not passing an LFMCMC object -----------------------
# Target is 56 tests
expected_error_msg <- "must be an object of class epiworld_lfmcmc"
not_lfmcmc <- c("NOT LFMCMC")

expect_error(run_lfmcmc(not_lfmcmc,
                        params_init = par0_int,
                        n_samples = n_samp_double,
                        epsilon = epsil_int,
                        seed = model_seed
                        ), expected_error_msg)

expect_error(set_observed_data(not_lfmcmc, obs_data), expected_error_msg)

expect_error(set_proposal_fun(not_lfmcmc, propfun), expected_error_msg)
expect_error(use_proposal_norm_reflective(not_lfmcmc), expected_error_msg)
expect_error(set_simulation_fun(not_lfmcmc, simfun), expected_error_msg)
expect_error(set_summary_fun(not_lfmcmc, sumfun), expected_error_msg)
expect_error(set_kernel_fun(not_lfmcmc, kernelfun), expected_error_msg)
expect_error(use_kernel_fun_gaussian(not_lfmcmc), expected_error_msg)

expect_error(set_params_names(not_lfmcmc, c("Par 1", "Par 2")), expected_error_msg)
expect_error(set_stats_names(not_lfmcmc, get_states(model_sir)), expected_error_msg)

expect_error(get_mean_params(not_lfmcmc), expected_error_msg)
expect_error(get_mean_stats(not_lfmcmc), expected_error_msg)
expect_error(get_accepted_params(not_lfmcmc), expected_error_msg)
expect_error(get_accepted_stats(not_lfmcmc), expected_error_msg)
expect_error(get_sample_stats(not_lfmcmc), expected_error_msg)
expect_error(get_n_params(not_lfmcmc), expected_error_msg)
expect_error(get_n_stats(not_lfmcmc), expected_error_msg)
expect_error(get_n_samples(not_lfmcmc), expected_error_msg)