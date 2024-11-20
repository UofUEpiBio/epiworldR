
set.seed(4222)
Y <- rnorm(2000, mean = -5, sd = 2.5)

simfun <- function(par) {
  rnorm(2000, mean = par[1], sd = par[2])
}

sumfun <- function(x) {
  c(mean(x), sd(x))
}

propfun <- function(par) {
  
  par_new <- par + rnorm(2, sd = 0.1)

  # Reflecting par2
  if (par_new[2] < 0) {
    par_new[2] <- par[2] - (par_new[2] - par[2])
  }

  return(par_new)
}

kernelfun <- function(stats_now, stats_obs, epsilon) {

  dnorm(sqrt(sum((stats_obs - stats_now)^2)))
  
}


lfmcmc_model <- LFMCMC() |>
  set_simulation_fun(simfun) |>
  set_summary_fun(sumfun) |>
  set_proposal_fun(propfun) |>
  set_kernel_fun(kernelfun) |>
  set_observed_data(Y)

x <- run_lfmcmc(
  lfmcmc = lfmcmc_model,
  params_init_ = c(0, 1),
  n_samples_ = 3000,
  epsilon_ = 1.0,
  seed = 4222
)

x_means <- get_params_mean(x)

expect_equivalent(
  x_means,
  c(-5, 2.5),
  tolerance = 0.1
)
