# Test just this file: tinytest::run_test_file("inst/tinytest/test-reproduction-number.R")

library(epiworldR)

mixing_matrix <- matrix(
  c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ),
  nrow = 3,
  byrow = TRUE
)

contact_rate <- 20
transmission_rate <- 0.2
recovery_rate <- 0.5
infectious_period_days <- 1 / recovery_rate
nsims <- 500

model <- ModelSEIRMixing(
  name = "R number check",
  n = 600,
  prevalence = 1 / 600,
  contact_rate = contact_rate,
  transmission_rate = transmission_rate,
  incubation_days = 3,
  recovery_rate = recovery_rate,
  contact_matrix = mixing_matrix
)

model |>
  add_entity(entity("Group 1", 200, FALSE)) |>
  add_entity(entity("Group 2", 200, FALSE)) |>
  add_entity(entity("Group 3", 200, FALSE))


verbose_off(model)

saver <- make_saver("reproductive")

run_multiple(
  model,
  ndays = 100,
  nsims = nsims,
  seed = 1234,
  saver = saver,
  verbose = FALSE,
  nthreads = 1L
)

entity_sizes <- vapply(get_entities(model), get_entity_size, numeric(1))
expect_equal(entity_sizes, rep(200, 3))

reproductive <- run_multiple_get_results(model, nthreads = 1L)$reproductive

# With prevalence = 1/600, the rows exposed on day 0 correspond to the single
# seeded case in each simulation, so their offspring count is the empirical Rt.
initial_rt <- reproductive$rt[
  reproductive$source != -1L &
    reproductive$source_exposure_date == 0L
]

expected_rt <- compute_reproduction_number(
  contact_matrix = contact_rate * mixing_matrix,
  transmission_prob = transmission_rate,
  infectious_period_days = infectious_period_days,
  infectiousness = c(1, 1, 1),
  susceptibility = c(1, 1, 1)
)

expect_equal(expected_rt$type, "R0")
expect_equal(length(initial_rt), nsims)
expect_true(abs(mean(initial_rt) - expected_rt$R) < 1.0)

# R0 = p_t * c / (1 - p_r)
R0_naive <- transmission_rate * contact_rate / recovery_rate
