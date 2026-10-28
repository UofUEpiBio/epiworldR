# Test just this file: tinytest::run_test_file("inst/tinytest/test-seir.R")

# Create small world population SEIR Model -------------------------------------
expect_silent(seir_0 <- ModelSEIR(
  name = "A Virus",
  n = 10000,
  prevalence = .01,
  transmission_prob = .5,
  latent_period = 7.0,
  recovery_rate = .1
))

# Check model initialization
expect_inherits(seir_0, "epiworld_seir")
expect_inherits(seir_0, "epiworld_model")
expect_length(class(seir_0), 2)
expect_silent(add_layer(seir_0, layer_network("network", graph_smallworld(k = 5, p = .01))))


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
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = bad_numeric_input,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = bad_numeric_input,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = bad_numeric_input,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = bad_numeric_input
), expected_error_msg_double)

# Check NA
expect_error(seir_0 <- ModelSEIR(
  name = NA,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = NA,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = NA,
  latent_period = good_incubation_days,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = NA,
  recovery_rate = good_recovery_rate
), expected_error_msg_na)

expect_error(seir_0 <- ModelSEIR(
  name = good_name,
  n          = 1000,
  prevalence = good_prevalence,
  transmission_prob = good_transmission_rate,
  latent_period = good_incubation_days,
  recovery_rate = NA
), expected_error_msg_na)