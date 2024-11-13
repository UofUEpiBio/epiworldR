# Test just this file: tinytest::run_test_file("inst/tinytest/test-sir.R")

# Create function to test transition probability matrix ------------------------
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
sir_0 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3
)

# Check model initialized correctly
expect_inherits(sir_0,"epiworld_sir")
expect_inherits(sir_0,"epiworld_model")
expect_silent(agents_smallworld(
  sir_0,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
))

# Running with queuing
# - Suppressing print output
expect_silent(verbose_off(sir_0))
# - Check plot fails before model is run
expect_warning(expect_error(plot(sir_0)))
# - Run model
expect_silent(run(sir_0, ndays = 50, seed = 1912))
# - Check plot succeeded after running model
expect_silent(plot(sir_0))

# Check transition probability matrix
tmat_0_queuing <- get_transition_probability(sir_0)
test_tmat_matches_expected(tmat_0_queuing)

# Run again without queuing and verify output is the same
expect_silent(queuing_off(sir_0))
run(sir_0, ndays = 50, seed = 1912)
tmat_0_noqueuing <- get_transition_probability(sir_0)
test_tmat_matches_expected(tmat_0_noqueuing)

expect_identical(tmat_0_noqueuing, tmat_0_queuing)

# Create new SIR model without queuing -----------------------------------------
sir_1 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3
)

agents_smallworld(
  sir_1,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
)

queuing_off(sir_1)

# Run the model and check output
verbose_off(sir_1)
run(sir_1, ndays = 50, seed = 1912)
tmat_1_noqueuing <- get_transition_probability(sir_1)

expect_identical(tmat_1_noqueuing, tmat_0_queuing)

# Check queuing is faster ------------------------------------------------------
runtime_0_noqueuing <- system.time(run(sir_0, ndays = 50, seed = 1912))
queuing_on(sir_0)
runtime_0_queuing <- system.time(run(sir_0, ndays = 50, seed = 1912))
expect_true(runtime_0_queuing["elapsed"] < runtime_0_noqueuing["elapsed"])
