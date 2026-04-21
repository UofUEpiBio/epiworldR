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
  transmission_rate = transmission_rate,
  incubation_days = 3,
  recovery_rate = recovery_rate,
  contact_matrix = mixing_matrix * contact_rate
)

model |>
  add_entity(entity("Group 1", 200, FALSE)) |>
  add_entity(entity("Group 2", 300, FALSE)) |>
  add_entity(entity("Group 3", 100, FALSE))

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
  susceptibility = c(1, 1, 1),
  group_sizes = entity_sizes
)

if (require(multigroup.vaccine)) {
  
  mgv_r_calc <- vaxrepnum(
    meaninf = infectious_period_days,
    popsize = entity_sizes,
    trmat = contact_rate * mixing_matrix * transmission_rate,
    initR = rep(0, 3),
    initV = rep(0, 3),
    vaxeff = 0
  )

  message("The calibrated R0 for multigroup.vaccine is: ", mgv_r_calc)

}

expect_equal(expected_rt$type, "R0")
expect_equal(length(initial_rt), nsims)
expect_true(abs(mean(initial_rt) - expected_rt$R) < 1.0)

# R0 = p_t * c / (1 - p_r)
R0_naive <- transmission_rate * contact_rate / recovery_rate
expect_equal(expected_rt$R, R0_naive)

weighted_contact_matrix <- matrix(
  c(
    10, 4, 1,
    2, 8, 3,
    2, 2, 6
  ),
  nrow = 3,
  byrow = TRUE
)
weighted_infectiousness <- c(1.0, 0.9, 1.1)
weighted_susceptibility <- c(0.8, 1.0, 0.7)
weighted_group_sizes <- c(100, 200, 50)

weighted_rt <- compute_reproduction_number(
  contact_matrix = weighted_contact_matrix,
  transmission_prob = 0.15,
  infectious_period_days = 4,
  infectiousness = weighted_infectiousness,
  susceptibility = weighted_susceptibility,
  group_sizes = weighted_group_sizes
)

expected_weighted_ngm <-
  0.15 * 4 *
    diag(weighted_infectiousness / weighted_group_sizes, nrow = 3, ncol = 3) %*%
      weighted_contact_matrix %*%
      diag(weighted_susceptibility * weighted_group_sizes, nrow = 3, ncol = 3)

expect_equal(weighted_rt$next_generation_matrix, expected_weighted_ngm)
expect_equal(weighted_rt$group_sizes, weighted_group_sizes)

balanced_by_size <- matrix(
  c(
    12, 4,
    2, 9
  ),
  nrow = 2,
  byrow = TRUE
)

expect_silent(
  compute_reproduction_number(
    contact_matrix = balanced_by_size,
    transmission_prob = 0.1,
    infectious_period_days = 3,
    group_sizes = c(100, 200),
    check_reciprocity = TRUE
  )
)

expect_warning(
  compute_reproduction_number(
    contact_matrix = balanced_by_size,
    transmission_prob = 0.1,
    infectious_period_days = 3,
    check_reciprocity = TRUE
  ),
  "reciprocity check"
)
