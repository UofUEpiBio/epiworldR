# Test just this file: tinytest::run_test_file("inst/tinytest/test-seir.R")

# Create small world population SEIR Model -------------------------------------
expect_silent(seir_0 <- ModelSEIR(
  name = "A Virus",
  prevalence = .01,
  transmission_rate = .5,
  incubation_days = 7.0,
  recovery_rate = .1
))

# Check model initialization
expect_inherits(seir_0, "epiworld_seir")
expect_inherits(seir_0, "epiworld_model")
expect_length(class(seir_0), 2)
expect_silent(agents_smallworld(
  seir_0,
  n = 10000,
  k = 5,
  d = FALSE,
  p = .01
))


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