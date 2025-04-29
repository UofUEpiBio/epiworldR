# Test just this file: tinytest::run_test_file("inst/tinytest/test-seir.R")

# Create small world population SEIR Model -------------------------------------
expect_silent(seir_0 <- ModelSEIR(
  name = "A Virus",
  prevalence = .01,
  transmission_rate = .5,
  incubation_days = 7.0,
  recovery_rate = .1
))

# Check model initialization
expect_inherits(seir_0, "epiworld_seir")
expect_inherits(seir_0, "epiworld_model")
expect_length(class(seir_0), 2)
expect_silent(agents_smallworld(
  seir_0,
  n = 10000,
  k = 5,
  d = FALSE,
  p = .01
))


# Check model run without queuing ----------------------------------------------
expect_silent(verbose_off(seir_0))
expect_silent(queuing_off(seir_0))
expect_silent(initial_states(seir_0, c(.3, .5)))
expect_error(plot(seir_0), "model must be run before it can be plotted")
expect_silent(run(seir_0, ndays = 100, seed = 1231))
expect_silent(plot(seir_0)) # Plot succeeds after model is run

hist_0 <- get_hist_total(seir_0)

expect_equal(hist_0[1,3], 4950)
expect_equal(hist_0[2,3], 70)
expect_equal(hist_0[3,3], 30)
expect_equal(hist_0[4,3], 4950)

# Check functions fail with invalid inputs -------------------------------------
good_name <- "A Virus"
good_prevalence <- .01
good_transmission_rate <- 0.5
good_incubation_days <- 4.0
good_recovery_rate <- 1.0/7.0

bad_name <- 10
bad_numeric_input <- "not a number"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(seir_0 <- ModelSEIR(
  name = bad_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = bad_numeric_input,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = bad_numeric_input,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = bad_numeric_input,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = bad_numeric_input
), expected_error_msg_double)

# Check NA
expect_error(seir_0 <- ModelSEIR(
  name = NA,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = NA,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = NA,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = NA,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = NA
), expected_error_msg_na)