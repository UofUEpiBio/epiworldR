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
expect_length(class(sird_0), 2)
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
expect_error(plot(sird_0), "model must be run before it can be plotted")
expect_silent(run(sird_0, ndays = 100, seed = 1231))
expect_silent(plot(sird_0)) # Plot succeeds after model is run

hist_0 <- get_hist_total(sird_0)

expect_equal(hist_0[1,3], 8931)
expect_equal(hist_0[2,3], 100)
expect_equal(hist_0[3,3], 474)
expect_equal(hist_0[4,3], 495)

# Check functions fail with invalid inputs -------------------------------------
good_name <- "A Virus"
good_prevalence <- .01
good_transmission_rate <- 0.9
good_recovery_rate <- 0.3
good_death_rate <- .1

bad_name <- NA
bad_prevalence <- NA
bad_transmission_rate <- NA
bad_recovery_rate <- NA
bad_death_rate <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(sird_0 <- ModelSIRD(
  name = bad_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_str)

expect_error(sird_0 <- ModelSIRD(
  name = good_name,
  prevalence = bad_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sird_0 <- ModelSIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = bad_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sird_0 <- ModelSIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = bad_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sird_0 <- ModelSIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = bad_death_rate
), expected_error_msg_double)