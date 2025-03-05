# Test just this file: tinytest::run_test_file("inst/tinytest/test-sir.R")

# Function to test transition probability matrix ------------------------
test_tmat_matches_expected <- function(tmat) {
  tmat_expected <- structure(
    c(
      0.961843, 0, 0,
      0.03815696, 0.69985167, 0,
      0, 0.3001483, 1
    ),
    dim = c(3L, 3L),
    dimnames = list(
      c("Susceptible", "Infected", "Recovered"),
      c("Susceptible", "Infected", "Recovered")
    )
  )
  
  expect_equal(tmat, tmat_expected, tolerance = 0.0000001)
}

# Create small world population SIR Model --------------------------------------
expect_silent(sir_0 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3
))

# Check model initialization
expect_inherits(sir_0, "epiworld_sir")
expect_inherits(sir_0, "epiworld_model")
expect_length(class(sir_0), 2)
expect_silent(agents_smallworld(
  sir_0,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
))

# Check model run with queuing -------------------------------------------------
expect_silent(verbose_off(sir_0))
expect_error(plot(sir_0), "model must be run before it can be plotted")
expect_silent(run(sir_0, ndays = 50, seed = 1912))
expect_silent(plot(sir_0)) # Plot succeeds after model is run

tmat_queuing <- get_transition_probability(sir_0)
test_tmat_matches_expected(tmat_queuing)

# Check model run without queuing ----------------------------------------------
expect_silent(queuing_off(sir_0))
run(sir_0, ndays = 50, seed = 1912)

tmat_noqueuing <- get_transition_probability(sir_0)
expect_identical(tmat_noqueuing, tmat_queuing)

# Check queuing is faster ------------------------------------------------------
runtime_noqueuing <- system.time(run(sir_0, ndays = 50, seed = 1912))
queuing_on(sir_0)
runtime_queuing <- system.time(run(sir_0, ndays = 50, seed = 1912))
expect_true(runtime_queuing["elapsed"] < runtime_noqueuing["elapsed"])

# Check functions fail with invalid inputs -------------------------------------
good_name <- "A Virus"
good_prevalence <- 0.01
good_transmission_rate <- 0.9
good_recovery_rate <- 0.3

bad_name <- NA
bad_prevalence <- NA
bad_transmission_rate <- NA
bad_recovery_rate <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(test_model <- ModelSIR(
  name = bad_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_str)

expect_error(test_model <- ModelSIR(
  name = good_name,
  prevalence = bad_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(test_model <- ModelSIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = bad_transmission_rate,
  recovery_rate = good_recovery_rate
), expected_error_msg_double)

expect_error(test_model <- ModelSIR(
  name = good_name,
  prevalence = good_prevalence,
  transmission_rate = good_transmission_rate,
  recovery_rate = bad_recovery_rate
), expected_error_msg_double)
