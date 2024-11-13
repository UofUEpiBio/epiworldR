# Test just this file: tinytest::run_test_file("inst/tinytest/test-sirconn.R")

# Function to test transition probability matrix ------------------------
test_tmat <- function(tmat) {
  tmat_expected <- structure(
    c(
      0.5397, 0, 0,
      0.4603, 0.8473519, 0,
      0, 0.1526481, 1
    ),
    dim = c(3L, 3L),
    dimnames = list(
      c("Susceptible", "Infected", "Recovered"),
      c("Susceptible", "Infected", "Recovered")
    )
  )
  
  # Check matches expected output
  expect_equal(tmat, tmat_expected, tolerance = 0.0000001)
  
  # Check for out of bounds values
  expect_false(any(tmat < 0))
  expect_false(any(tmat > 1))
}

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

# Check model run with queuing -------------------------------------------------
expect_silent(verbose_off(sirconn_0))
expect_warning(expect_error(plot(sirconn_0))) # Plot fails before model is run
expect_silent(run(sirconn_0, ndays = 100, seed = 131))
expect_silent(plot(sirconn_0)) # Plot succeeds after model is run

hist_queuing <- get_hist_total(sirconn_0)
tmat_queuing <- get_transition_probability(sirconn_0)

test_tmat(tmat_queuing)

# Check model run without queuing ----------------------------------------------
expect_silent(queuing_off(sirconn_0))
run(sirconn_0, ndays = 100, seed = 131)

hist_noqueuing <- get_hist_total(sirconn_0)
tmat_noqueuing <- get_transition_probability(sirconn_0)

expect_identical(hist_noqueuing, hist_queuing)
expect_identical(tmat_noqueuing, tmat_queuing)

