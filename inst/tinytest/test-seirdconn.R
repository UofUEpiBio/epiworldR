# Test just this file: tinytest::run_test_file("inst/tinytest/test-seirdconn.R")

# Create SEIRDCONN Model -------------------------------------------------------

good_name <- "A Virus"
good_n <- 10000
good_prevalence <- .01
good_contact_rate <- 2
good_incubation_days <- 7
good_transmission_rate <- 0.5
good_recovery_rate <- 0.3
good_death_rate <- 0.01

expect_silent(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
))

# Check model initialization
expect_inherits(seirdconn_0, "epiworld_seirdconn")
expect_inherits(seirdconn_0, "epiworld_model")
expect_length(class(seirdconn_0), 2)


# Check functions fail with invalid inputs -------------------------------------
bad_name <- 10
bad_numeric_input <- "not a number"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = bad_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_str)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = bad_numeric_input,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_int)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = bad_numeric_input,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = bad_numeric_input,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = bad_numeric_input,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = bad_numeric_input,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = bad_numeric_input,
  death_rate = good_death_rate
), expected_error_msg_double)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = bad_numeric_input
), expected_error_msg_double)

# Check NA

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = NA,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_str)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = NA,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = NA,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = NA,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = NA,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = NA,
  recovery_rate = good_recovery_rate,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = NA,
  death_rate = good_death_rate
), expected_error_msg_na)

expect_error(seirdconn_0 <- ModelSEIRDCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  incubation_days = good_incubation_days,
  recovery_rate = good_recovery_rate,
  death_rate = NA
), expected_error_msg_na)