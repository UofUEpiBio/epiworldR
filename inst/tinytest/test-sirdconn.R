# Test just this file: tinytest::run_test_file("inst/tinytest/test-sirdconn.R")

# Create SIRDCONN Model --------------------------------------------------------
good_name <- "A Virus"
good_n <- 10000
good_prevalence <- .01
good_contact_rate <- 5
good_transmission_rate <- 0.4
good_recovery_rate <- 0.5
good_death_rate <- 0.01

expect_silent(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
))

# Check model initialization
expect_inherits(sirdconn_0, "epiworld_sirdconn")
expect_inherits(sirdconn_0, "epiworld_model")
expect_length(class(sirdconn_0), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name <- 10
bad_numeric_input <- "not a number"

expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = bad_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_str)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = bad_numeric_input,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_int)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = bad_numeric_input,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = bad_numeric_input,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = bad_numeric_input,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = bad_numeric_input,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(sirdconn_0 <- ModelSIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate,
  death_rate = bad_numeric_input
), expected_error_msg_double)