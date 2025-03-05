# Test just this file: tinytest::run_test_file("inst/tinytest/test-sisd.R")

# Create SISD Model ------------------------------------------------------------
good_name               <- "A Virus"
good_prevalence         <- 0.01
good_transmission_rate  <- 0.9
good_recovery_rate      <- 0.1
good_death_rate         <- 0.01

expect_silent(sisd_0 <- ModelSISD(
  name                  = good_name,
  prevalence            = good_prevalence,
  transmission_rate     = good_transmission_rate,
  recovery_rate         = good_recovery_rate,
  death_rate            = good_death_rate
))

# Check model initialization
expect_inherits(sisd_0, "epiworld_sisd")
expect_inherits(sisd_0, "epiworld_model")
expect_length(class(sisd_0), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name                <- NA
bad_prevalence          <- NA
bad_transmission_rate   <- NA
bad_recovery_rate       <- NA
bad_death_rate          <- NA

expected_error_msg_str    <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(sisd_0 <- ModelSISD(
  name                = bad_name,
  prevalence          = good_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = good_recovery_rate,
  death_rate          = good_death_rate
), expected_error_msg_str)

expect_error(sisd_0 <- ModelSISD(
  name                = good_name,
  prevalence          = bad_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = good_recovery_rate,
  death_rate          = good_death_rate
), expected_error_msg_double)

expect_error(sisd_0 <- ModelSISD(
  name                = good_name,
  prevalence          = good_prevalence,
  transmission_rate   = bad_transmission_rate,
  recovery_rate       = good_recovery_rate,
  death_rate          = good_death_rate
), expected_error_msg_double)

expect_error(sisd_0 <- ModelSISD(
  name                = good_name,
  prevalence          = good_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = bad_recovery_rate,
  death_rate          = good_death_rate
), expected_error_msg_double)

expect_error(sisd_0 <- ModelSISD(
  name                = good_name,
  prevalence          = good_prevalence,
  transmission_rate   = good_transmission_rate,
  recovery_rate       = good_recovery_rate,
  death_rate          = bad_death_rate
), expected_error_msg_double)

