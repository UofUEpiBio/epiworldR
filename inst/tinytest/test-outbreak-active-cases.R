# Test just this file: tinytest::run_test_file("inst/tinytest/test-outbreak-active-cases.R")

# Test get_outbreak_size and get_active_cases functions
# These tests verify that:
# 1. get_active_cases() returns the number of infected individuals
# 2. get_outbreak_size() returns the number of infected + recovered individuals

# Create an SIR model with random mixing ---------------------------------------
model <- ModelSIR(
  name = "test virus",
  n = 100,
  prevalence = 0.1,
  transmission_prob = 0.5,
  recovery_rate = 0.3
)
add_param(model, "Contact rate", 2)
add_layer(
  model,
  layer_mixing("everyone", groups = rep(0, 100), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

verbose_off(model)
run(model, ndays = 15, seed = 123)

# Get model history data (for comparison)
hist_total <- get_hist_total(model)

# Get active cases and outbreak size
active_cases <- get_active_cases(model)
outbreak_size <- get_outbreak_size(model)

# Checking active cases
hist_total_infected <- aggregate(
  counts ~ date,
  data = hist_total[hist_total$state == "Infected", ],
  sum
)

expect_equal(
  hist_total_infected$counts,
  active_cases$active_cases
)

# Checking outbreak size
hist_total_outbreak_size <- aggregate(
  counts ~ date,
  data = hist_total[hist_total$state %in% c("Infected", "Recovered"), ],
  sum
)

expect_equal(
  hist_total_outbreak_size$counts,
  outbreak_size$outbreak_size
)
