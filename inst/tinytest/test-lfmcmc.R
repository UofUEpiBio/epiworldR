# Create Model to use in LFMCMC simulation
model_seed <- 122

model_sir <- ModelSIR(
  name = "COVID-19",
  prevalence = .1,
  transmission_rate = .9,
  recovery_rate = .3
)

agents_smallworld(
  model_sir,
  n = 1000,
  k = 5,
  d = FALSE,
  p = 0.01
)

verbose_off(model_sir)

run(
  model_sir,
  ndays = 50,
  seed = model_seed
)

# Setup LFMCMC
## Extract the observed data from the model
obs_data <- unname(as.integer(get_today_total(model_sir)))

## Define the LFMCMC simulation function
simfun <- function(params) {

  set_param(model_sir, "Recovery rate", params[1])
  set_param(model_sir, "Transmission rate", params[2])

  run(
    model_sir,
    ndays = 50
  )

  res <- unname(as.integer(get_today_total(model_sir)))
  return(res)
}

## Define the LFMCMC summary function
sumfun <- function(dat) {
  return(dat)
}

## Define the LFMCMC proposal function
## - Based on proposal_fun_normal from lfmcmc-meat.hpp
propfun <- function(params_prev) {
  res <- params_prev + rnorm(length(params_prev), )
  return(res)
}

## Define the LFMCMC kernel function
## - Based on kernel_fun_uniform from lfmcmc-meat.hpp
kernelfun <- function(stats_now, stats_obs, epsilon) {

  ans <- sum(mapply(function(v1, v2) (v1 - v2)^2,
    stats_obs,
    stats_now))

  return(ifelse(sqrt(ans) < epsilon, 1.0, 0.0))
}

## Create the LFMCMC model
lfmcmc_model <- LFMCMC(model_sir) |>
  set_simulation_fun(simfun) |>
  set_summary_fun(sumfun) |>
  set_proposal_fun(propfun) |>
  set_kernel_fun(kernelfun) |>
  set_observed_data(obs_data)

# Run LFMCMC simulation
## Set initial parameters
par0 <- as.double(c(0.1, 0.5))
n_samp <- 2000
epsil <- as.double(1.0)

## Run the LFMCMC simulation
run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = par0,
  n_samples_ = n_samp,
  epsilon_ = epsil,
  seed = model_seed
)

# Print the results
set_stats_names(lfmcmc_model, get_states(model_sir))
set_par_names(lfmcmc_model, c("Immune recovery", "Infectiousness"))

print(lfmcmc_model)