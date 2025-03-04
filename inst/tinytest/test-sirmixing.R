# Test just this file: tinytest::run_test_file("inst/tinytest/test-sirmixing.R")

good_name               <- "A Virus"
good_n                  <- 9e3
good_prevalence         <- 1 / good_n
good_contact_rate       <- 20
good_transmission_rate  <- 0.0
good_recovery_rate      <- 1 / 7
good_contact_matrix     <- c(
  c(0.9, 0.05, 0.05),
  c(0.1, 0.8, 0.1),
  c(0.1, 0.2, 0.7)
) |> matrix(byrow = TRUE, nrow = 3)

expect_silent(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
))

# Check model initialization
expect_inherits(test_model, "epiworld_sirmixing")
expect_inherits(test_model, "epiworld_model")
expect_length(class(test_model), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name               <- NA
bad_n                  <- NA
bad_prevalence         <- NA
bad_contact_rate       <- NA
bad_transmission_rate  <- NA
bad_recovery_rate      <- NA
bad_contact_matrix     <- NA

expected_error_msg_str    <- "must be a string"
expected_error_msg_int    <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

expect_error(test_model <- ModelSIRMixing(
  name              = bad_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_str)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = bad_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_int)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = bad_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = bad_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = bad_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = bad_recovery_rate,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  contact_matrix    = bad_contact_matrix
), expected_error_msg_any_na)