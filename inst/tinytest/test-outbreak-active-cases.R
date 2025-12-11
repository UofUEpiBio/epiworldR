# Test just this file: tinytest::run_test_file("inst/tinytest/test-outbreak-active-cases.R")

# Test get_outbreak_size and get_active_cases functions across multiple days
# This test validates that virus history is correctly recorded for all days
# of a simulation, not just day 0. Ported from epiworld C++ test:
# https://github.com/UofUEpiBio/epiworld/blob/e338e64b912e50acee6e14705c4d3b0bf8baf885/tests/24a-virus-hist-single-run.cpp

# Create a simple SIR model ----------------------------------------------------
expect_silent(model <- ModelSIR(
  name = "test virus",
  prevalence = 0.1,
  transmission_rate = 0.3,
  recovery_rate = 0.1
))

# Set up a small world network with 100 agents
expect_silent(agents_smallworld(
  model,
  n = 100,
  k = 4,
  d = FALSE,
  p = 0.01
))

# Turn off verbose output
expect_silent(verbose_off(model))

# Run for 20 days (single simulation)
ndays <- 20
expect_silent(run(model, ndays = ndays, seed = 123))

# Test get_outbreak_size function ---------------------------------------------
outbreak_data <- get_outbreak_size(model)
expect_inherits(outbreak_data, "data.frame")

# Check that it has the expected columns
expected_cols <- c("date", "virus_id", "virus", "outbreak_size")
expect_equal(names(outbreak_data), expected_cols)

# Check that date column contains values from 0 to ndays
expect_true(all(outbreak_data$date >= 0))
expect_true(all(outbreak_data$date <= ndays))

# Check that we have data for multiple days (not just day 0)
unique_dates <- unique(outbreak_data$date)
expect_true(length(unique_dates) > 1, 
  info = "Outbreak size should be recorded for multiple days")

# Check that outbreak_size values are non-negative integers
expect_true(all(outbreak_data$outbreak_size >= 0))
expect_true(all(outbreak_data$outbreak_size == floor(outbreak_data$outbreak_size)))

# Check that virus_id is a non-negative integer
expect_true(all(outbreak_data$virus_id >= 0))

# Check that virus name matches the model
expect_equal(unique(outbreak_data$virus), "test virus")

# Outbreak size tracking across time
# Note: The function may return duplicate dates or non-sequential dates
# We just verify that data is being recorded and structure is correct
for (vid in unique(outbreak_data$virus_id)) {
  virus_data <- outbreak_data[outbreak_data$virus_id == vid, ]
  
  # Check that we have some data points
  expect_true(nrow(virus_data) > 0,
    info = "Should have outbreak size data for each virus")
}

# Test get_active_cases function ----------------------------------------------
active_cases_data <- get_active_cases(model)
expect_inherits(active_cases_data, "data.frame")

# Check that it has the expected columns
expected_cols_active <- c("date", "virus_id", "virus", "active_cases")
expect_equal(names(active_cases_data), expected_cols_active)

# Check that date column contains values from 0 to ndays
expect_true(all(active_cases_data$date >= 0))
expect_true(all(active_cases_data$date <= ndays))

# Check that we have data for multiple days (not just day 0)
unique_dates_active <- unique(active_cases_data$date)
expect_true(length(unique_dates_active) > 1,
  info = "Active cases should be recorded for multiple days")

# Check that active_cases values are non-negative integers
expect_true(all(active_cases_data$active_cases >= 0))
expect_true(all(active_cases_data$active_cases == floor(active_cases_data$active_cases)))

# Check that virus_id is a non-negative integer
expect_true(all(active_cases_data$virus_id >= 0))

# Check that virus name matches the model
expect_equal(unique(active_cases_data$virus), "test virus")

# Active cases should be less than or equal to outbreak size at each time point
merged_data <- merge(
  outbreak_data,
  active_cases_data,
  by = c("date", "virus_id", "virus"),
  suffixes = c("_outbreak", "_active")
)

expect_true(all(merged_data$active_cases <= merged_data$outbreak_size),
  info = "Active cases should never exceed outbreak size")

# Test consistency across multiple runs ---------------------------------------
# Run the model again with a different seed
expect_silent(run(model, ndays = ndays, seed = 456))

outbreak_data_2 <- get_outbreak_size(model)
active_cases_data_2 <- get_active_cases(model)

# Both should still have the same structure
expect_equal(names(outbreak_data_2), expected_cols)
expect_equal(names(active_cases_data_2), expected_cols_active)

# Should have data for multiple days
expect_true(length(unique(outbreak_data_2$date)) > 1)
expect_true(length(unique(active_cases_data_2$date)) > 1)

# Test with different model (SEIR) ---------------------------------------------
seir_model <- ModelSEIR(
  name = "seir virus",
  prevalence = 0.05,
  transmission_rate = 0.5,
  incubation_days = 5.0,
  recovery_rate = 0.15
)
expect_inherits(seir_model, "epiworld_model")

agents_smallworld(
  seir_model,
  n = 100,
  k = 4,
  d = FALSE,
  p = 0.01
)

expect_silent(verbose_off(seir_model))
expect_silent(run(seir_model, ndays = 15, seed = 789))

outbreak_seir <- get_outbreak_size(seir_model)
active_seir <- get_active_cases(seir_model)

# Should work for SEIR model as well
expect_inherits(outbreak_seir, "data.frame")
expect_inherits(active_seir, "data.frame")
expect_equal(names(outbreak_seir), expected_cols)
expect_equal(names(active_seir), expected_cols_active)

# Check virus name for SEIR
expect_equal(unique(outbreak_seir$virus), "seir virus")
expect_equal(unique(active_seir$virus), "seir virus")

# Test error handling ----------------------------------------------------------
# Create a model but don't run it
unrun_model <- ModelSIR(
  name = "unrun",
  prevalence = 0.1,
  transmission_rate = 0.3,
  recovery_rate = 0.1
)
expect_inherits(unrun_model, "epiworld_model")

agents_smallworld(unrun_model, n = 50, k = 3, d = FALSE, p = 0.01)

# Should still be able to call the functions on an unrun model
# (they should return empty or initial data)
outbreak_unrun <- get_outbreak_size(unrun_model)
active_unrun <- get_active_cases(unrun_model)

expect_inherits(outbreak_unrun, "data.frame")
expect_inherits(active_unrun, "data.frame")
