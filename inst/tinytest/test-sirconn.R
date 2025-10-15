# Test just this file: tinytest::run_test_file("inst/tinytest/test-sirconn.R")

# Create SIR CONN Model --------------------------------------------------------
expect_silent(sirconn_0 <- ModelSIRCONN(
  name = "A Virus",
  n = 10000,
  prevalence = .01,
  contact_rate = 4.0,
  transmission_rate = .5,
  recovery_rate = 1.0/7.0
))

# Check model initialization
expect_inherits(sirconn_0, "epiworld_sirconn")
expect_inherits(sirconn_0, "epiworld_model")
expect_length(class(sirconn_0), 2)

# Check model run with queuing -------------------------------------------------
expect_silent(verbose_off(sirconn_0))
expect_warning(queuing_on(sirconn_0), "SIR Connected models do not have queue.")
expect_error(plot(sirconn_0), "model must be run before it can be plotted")
expect_silent(run(sirconn_0, ndays = 100, seed = 131))
expect_silent(plot(sirconn_0)) # Plot succeeds after model is run

hist_queuing <- get_hist_total(sirconn_0)
tmat_queuing <- get_transition_probability(sirconn_0)

test_tmat_matches_expected(tmat_queuing)

# Check model run without queuing ----------------------------------------------
expect_warning(queuing_off(sirconn_0), "SIR Connected models do not have queue.")
run(sirconn_0, ndays = 100, seed = 131)

hist_noqueuing <- get_hist_total(sirconn_0)
tmat_noqueuing <- get_transition_probability(sirconn_0)

expect_identical(hist_noqueuing, hist_queuing)
expect_identical(tmat_noqueuing, tmat_queuing)

# Check functions fail with invalid inputs -------------------------------------
good_name <- "A Virus"
good_n <- 10000
good_prevalence <- .01
good_contact_rate <- 4.0
good_transmission_rate <- 0.5
good_recovery_rate <- 1.0/7.0

bad_name <- 10
bad_numeric_input <- "not a number"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

expect_error(sirconn_0 <- ModelSIRCONN(
  name = bad_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(sirconn_0 <- ModelSIRCONN(
  name = good_name,
  n = bad_numeric_input,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_int)

expect_error(sirconn_0 <- ModelSIRCONN(
  name = good_name,
  n = good_n,
  prevalence = bad_numeric_input,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(sirconn_0 <- ModelSIRCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = bad_numeric_input,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(sirconn_0 <- ModelSIRCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = bad_numeric_input,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(sirconn_0 <- ModelSIRCONN(
  name = good_name,
  n = good_n,
  prevalence = good_prevalence,
  contact_rate = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate = bad_numeric_input
), expected_error_msg_double)

