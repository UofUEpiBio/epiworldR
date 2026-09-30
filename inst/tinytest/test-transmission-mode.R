# Test just this file: tinytest::run_test_file("inst/tinytest/test-transmission-mode.R")

make_model <- function() {
  model <- ModelSIR(
    name = "A Virus", prevalence = .01, transmission_rate = .3,
    recovery_rate = 1/7
  )
  agents_smallworld(model, n = 5000, k = 5, d = FALSE, p = .01)
  verbose_off(model)
  model
}

# Defaults ---------------------------------------------------------------------
model <- make_model()
expect_equal(get_transmission_mode(model), "auto")
expect_equal(get_transmission_kappa(model), 0.5)

# Forced modes are used on every day, and a seed reproduces a run ------------
hist_of <- function(mode) {
  model <- make_model()
  expect_identical(set_transmission_mode(model, mode), model)
  expect_equal(get_transmission_mode(model), mode)
  run(model, ndays = 40, seed = 1912)
  expect_equal(get_last_transmission_mode(model), mode)
  h <- get_hist_total(model)
  h[order(h$date, h$state), ]
}

push_1 <- hist_of("push")
push_2 <- hist_of("push")
pull_1 <- hist_of("pull")

expect_identical(push_1, push_2)

# Same distribution, different random numbers: the runs differ, but both
# produce an outbreak.
expect_false(identical(push_1, pull_1))
final_recovered <- function(h) h$counts[h$date == 40 & h$state == "Recovered"]
expect_true(final_recovered(push_1) > 50)
expect_true(final_recovered(pull_1) > 50)

# kappa ------------------------------------------------------------------------
model <- make_model()
set_transmission_mode(model, "auto", kappa = 2)
expect_equal(get_transmission_kappa(model), 2)

# kappa = 0 means "auto" pulls whenever there are carriers to push from
set_transmission_mode(model, "auto", kappa = 0)
run(model, ndays = 10, seed = 1)
expect_equal(get_last_transmission_mode(model), "pull")

# NULL goes back to the library default
set_transmission_mode(model, "auto")
expect_equal(get_transmission_kappa(model), 0.5)

# Models without the default network sampler always pull ----------------------
mixing <- ModelSIRMixing(
  name = "Flu", n = 1000, prevalence = .01, transmission_rate = .1,
  recovery_rate = 1/7, contact_matrix = matrix(5, 1, 1)
)
mixing |> add_entity(entity("All", 1000, as_proportion = FALSE))
verbose_off(mixing)
set_transmission_mode(mixing, "push")
run(mixing, ndays = 5, seed = 1)
expect_equal(get_last_transmission_mode(mixing), "pull")

# Invalid inputs ---------------------------------------------------------------
model <- make_model()
expect_error(set_transmission_mode(model, "sideways"), "should be one of")
expect_error(set_transmission_mode(model, "push", kappa = -1), "non-negative")
expect_error(set_transmission_mode(model, "push", kappa = NA), "non-negative")
expect_error(set_transmission_mode(model, "push", kappa = c(1, 2)), "non-negative")
expect_error(set_transmission_mode(NULL, "push"), "epiworld_model")
expect_error(get_transmission_mode(NULL), "epiworld_model")
