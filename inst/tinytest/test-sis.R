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
expect_length(class(sis_0), 2)
expect_silent(agents_smallworld(
  sis_0,
  n = 2000,
  k = 5,
  d = FALSE,
  p = .01
))

# Check model run --------------------------------------------------------------
expect_silent(verbose_off(sis_0))
expect_error(plot(sis_0), "model must be run before it can be plotted")
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

# Check functions fail with invalid inputs -------------------------------------
good_name               <- "A Virus"
good_prevalence         <- 0.1
good_transmission_rate  <- 0.3
good_recovery_rate      <- 0.3

bad_name                <- NA
bad_prevalence          <- NA
bad_transmission_rate   <- NA
bad_recovery_rate       <- NA

expected_error_msg_str    <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(sis_0 <- ModelSIS(
  name                = bad_name,
  prevalence          = good_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = good_recovery_rate
), expected_error_msg_str)

expect_error(sis_0 <- ModelSIS(
  name                = good_name,
  prevalence          = bad_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = good_recovery_rate
), expected_error_msg_double)

expect_error(sis_0 <- ModelSIS(
  name                = good_name,
  prevalence          = good_prevalence,
  transmission_rate   = bad_transmission_rate,
  recovery_rate       = good_recovery_rate
), expected_error_msg_double)

expect_error(sis_0 <- ModelSIS(
  name                = good_name,
  prevalence          = good_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = bad_recovery_rate
), expected_error_msg_double)
