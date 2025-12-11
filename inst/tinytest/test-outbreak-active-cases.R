# Test just this file: tinytest::run_test_file("inst/tinytest/test-outbreak-active-cases.R")

# Test get_outbreak_size and get_active_cases functions
# These tests verify that:
# 1. get_active_cases() returns the number of infected individuals
# 2. get_outbreak_size() returns the number of infected + recovered individuals

# Create a SIRCONN model -------------------------------------------------------
model <- ModelSIRCONN(
  name = "test virus",
  n = 100,
  prevalence = 0.1,
  contact_rate = 2,
  transmission_rate = 0.5,
  recovery_rate = 0.3
)

verbose_off(model)
run(model, ndays = 15, seed = 123)

# Get model history data (for comparison)
hist_total <- get_hist_total(model)

# Get active cases and outbreak size
active_cases <- get_active_cases(model)
outbreak_size <- get_outbreak_size(model)

# Test 1: get_active_cases() should match number of infected -------------------
# For each day, active cases should equal the infected count from history

for (day in unique(active_cases$date)) {
  # Get infected count from history
  infected_count <- hist_total[hist_total$date == day & hist_total$state == "Infected", "counts"]
  
  # Get active cases for this day (should be only one row per day)
  active_count <- active_cases[active_cases$date == day, "active_cases"]
  
  # They should be equal
  expect_equal(active_count, infected_count,
    info = paste("Active cases should equal infected count on day", day))
}

# Test 2: get_outbreak_size() should match infected + recovered ---------------
# For each day in outbreak_size, outbreak size should equal infected + recovered 
# from history. Note: we deduplicate by taking the first occurrence of each date.

outbreak_size_unique <- outbreak_size[!duplicated(outbreak_size$date), ]

for (i in seq_len(nrow(outbreak_size_unique))) {
  day <- outbreak_size_unique$date[i]
  outbreak_count <- outbreak_size_unique$outbreak_size[i]
  
  # Get infected and recovered counts from history
  infected_count <- hist_total[hist_total$date == day & hist_total$state == "Infected", "counts"]
  recovered_count <- hist_total[hist_total$date == day & hist_total$state == "Recovered", "counts"]
  
  # Outbreak size should equal infected + recovered
  expect_equal(outbreak_count, infected_count + recovered_count,
    info = paste("Outbreak size should equal infected + recovered on day", day))
}
