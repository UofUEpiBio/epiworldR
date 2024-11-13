# Run just this file with command:
# - tinytest::run_test_file("inst/tinytest/test-sir.R")
# Adding a Small world population without queuing ------------------------------
sir_0 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3
)

# Check model initialized correctly
expect_inherits(sir_0,"epiworld_sir")
expect_inherits(sir_0,"epiworld_model")
# Check can't plot before model is run
expect_warning(expect_error(plot(sir_0)))

expect_silent(agents_smallworld(
  sir_0,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
))

# Initializing 
queuing_off(sir_0)

# Running
verbose_off(sir_0)
run(sir_0, ndays = 50, seed = 1912)

# Check plots without issue after running model
expect_silent(plot(sir_0))

tmat_0 <- get_transition_probability(sir_0)

expected_dimnames <- list(
  c("Susceptible", "Infected", "Recovered"),
  c("Susceptible", "Infected", "Recovered")
)
expected_dim <- c(3L, 3L)

expect_equal(dimnames(tmat_0), expected_dimnames)
expect_equal(dim(tmat_0), expected_dim)

# Creating a SIR model with queuing --------------------------------------------
# TODO: Can we simply reset the model and run again?
sir_1 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  transmission_rate = .9,
  recovery_rate = .3
)

# Adding a Small world population 
agents_smallworld(
  sir_1,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
)

# Running and printing
verbose_off(sir_1)
run(sir_1, ndays = 50, seed = 1912)

tmat_1 <- get_transition_probability(sir_1)

# Expected
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

# Check matches expected output
expect_equal(tmat_0, tmat_expected, tolerance = 0.0000001)
# Check SIR produces same output with and without queuing
# TODO: why is this the case?
expect_equal(tmat_0, tmat_1)