# Test just this file: tinytest::run_test_file("inst/tinytest/test-seird.R")

# Create SEIRD Model -----------------------------------------------------------

good_name <- "A Virus"
good_prevalence <- .01
good_transmission_rate <- 0.5
good_incubation_days <- 4.0
good_recovery_rate <- 1.0/7.0
good_death_rate <- 0.01

expect_silent(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
))

# Check model initialization
expect_inherits(seird_0, "epiworld_seird")
expect_inherits(seird_0, "epiworld_model")
expect_length(class(seird_0), 2)


# Check functions fail with invalid inputs -------------------------------------
bad_name <- NA
bad_prevalence <- NA
bad_transmission_rate <- NA
bad_incubation_days <- NA
bad_recovery_rate <- NA
bad_death_rate <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(seird_0 <- ModelSEIRD(
  name = bad_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_str)

expect_error(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = bad_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = bad_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = bad_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = bad_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seird_0 <- ModelSEIRD(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = bad_death_rate
), expected_error_msg_double)