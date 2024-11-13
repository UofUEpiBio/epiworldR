# Test just this file: tinytest::run_test_file("inst/tinytest/test-sird.R")

# Create SIRD Model ------------------------------------------------------------
expect_silent(sird_0 <- ModelSIRD(
  name = "A Virus",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3,
  death_rate = .1
))

# Check model initialization
expect_inherits(sird_0, "epiworld_sird")
expect_inherits(sird_0, "epiworld_model")
expect_silent(agents_smallworld(
  sird_0,
  n = 10000,
  k = 5,
  d = FALSE,
  p = .01
))

# Check model run --------------------------------------------------------------
expect_silent(verbose_off(sird_0))
expect_silent(initial_states(sird_0, c(.05, .05)))
expect_warning(expect_error(plot(sird_0))) # Plot fails before model is run
expect_silent(run(sird_0, ndays = 100, seed = 1231))
expect_silent(plot(sird_0)) # Plot succeeds after model is run

hist_0 <- get_hist_total(sird_0)

expect_equal(hist_0[1,3], 8931)
expect_equal(hist_0[2,3], 100)
expect_equal(hist_0[3,3], 474)
expect_equal(hist_0[4,3], 495)