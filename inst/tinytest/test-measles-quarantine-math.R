
# An in a school with low vaccination
R0s <- c(0.8, 1.25, 3.5)
for (R0 in R0s) {
  p_t <- .05
  p_r <- 1/7

  crate <- R0 / p_t * p_r

  model_measles <- ModelMeaslesQuarantine(
    n = 1000,
    prevalence = 1,
    contact_rate = crate,
    transmission_rate = p_t,
    prodromal_period = 1/p_r,
    prop_vaccinated = 0,
    quarantine_period = -1
  ) 

  # Running and printing
  saver <- make_saver("reproductive")
  run_multiple(
    model_measles, ndays = 60,
    seed = 1912,
    saver = saver,
    nsims = 200,
    nthreads = 2L
    )

  res <- run_multiple_get_results(model_measles, nthreads = 2L)

  # Identifying the date 0
  r0s <- res$reproductive
  r0s <- subset(r0s, source_exposure_date == 0 & source != -1)
  r0_obs <- r0s$rt |> mean()

  print(expect_equal(r0_obs, R0, tolerance = 0.1))

}