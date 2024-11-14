# Test just this file: tinytest::run_test_file("inst/tinytest/test-sis.R")

# Create small world population SIS Model --------------------------------------
expect_silent(sis_0 <- ModelSIS(
  name = "SIS",
  prevalence = .1,
  transmission_rate = .3,
  recovery_rate = .3
))

# Check model initialization
expect_inherits(sis_0, "epiworld_sis")
expect_inherits(sis_0, "epiworld_model")
expect_silent(agents_smallworld(
  sis_0,
  n = 2000,
  k = 5,
  d = FALSE,
  p = .01
))

# Check model run --------------------------------------------------------------
expect_silent(verbose_off(sis_0))
expect_warning(expect_error(plot(sis_0))) # Plot fails before model is run
expect_silent(run_multiple(
  sis_0,
  ndays=100,
  nsims=10,
  seed=1231,
  reset=TRUE,
  verbose = TRUE,
  nthreads = 1
))
expect_silent(plot(sis_0)) # Plot succeeds after model is run
