# Test just this file: tinytest::run_test_file("inst/tinytest/test-virus.R")

# Check creating valid distribution function
expect_silent(distfun <- distribute_virus_randomly(
  prevalence = 2L,
  as_proportion = FALSE,
  agents_ids = integer(0)
))

expect_inherits(distfun, "epiworld_virus_distfun")
expect_length(class(distfun), 1)

# Check distributing virus to invalid agent IDs
expect_error(distribute_virus_randomly(
  prevalence = 2L,
  as_proportion = FALSE,
  agents_ids = c(-1, -2, -33)
), "ID must be a positive integer")

# Test virus creation and printing ---------------------------------------------
expect_silent(virus_1 <- virus(
  name = "Test Virus",
  prevalence = 0.01,
  as_proportion = TRUE,
  prob_infecting = 0.5,
  recovery_rate = 0.3
))
expect_inherits(virus_1, "epiworld_virus")
expect_stdout(print(virus_1))

# Test virus probability setters -----------------------------------------------
expect_silent(set_prob_infecting(virus_1, 0.8))
expect_silent(set_prob_recovery(virus_1, 0.4))
expect_silent(set_prob_death(virus_1, 0.01))
expect_silent(set_incubation(virus_1, 5.0))

# Test virus name functions ----------------------------------------------------
expect_equal(get_name_virus(virus_1), "Test Virus")
expect_silent(set_name_virus(virus_1, "Modified Virus"))
expect_equal(get_name_virus(virus_1), "Modified Virus")

# Scenario (a): High infectiousness, no recovery ------------------------------
# With perfect transmission and no recovery, epidemic should spread widely
expect_silent(model_a <- ModelSIRCONN(
  name = "Virus Test (a)",
  n = 100,
  prevalence = 0.01,
  contact_rate = 10,
  transmission_rate = 1.0,
  recovery_rate = 0.0
))

expect_silent(verbose_off(model_a))
expect_silent(run(model_a, ndays = 10, seed = 123))

hist_a <- get_hist_total(model_a)
final_data_a <- hist_a[hist_a$date == max(hist_a$date), ]
final_infected_a <- final_data_a[final_data_a$state == "Infected", "counts"]
final_recovered_a <- final_data_a[final_data_a$state == "Recovered", "counts"]

expect_true(final_infected_a > 50, info = "Most agents should be infected")
expect_equal(final_recovered_a, 0, info = "No one should recover with 0 recovery rate")

# Scenario (b): No infectiousness, immediate recovery -------------------------
# Only the initially infected person should recover immediately
expect_silent(model_b <- ModelSIRCONN(
  name = "Virus Test (b)",
  n = 100,
  prevalence = 0.01,
  contact_rate = 10,
  transmission_rate = 0.0,
  recovery_rate = 1.0
))

expect_silent(verbose_off(model_b))
expect_silent(run(model_b, ndays = 1, seed = 123))

hist_b <- get_hist_total(model_b)
final_data_b <- hist_b[hist_b$date == max(hist_b$date), ]
final_infected_b <- final_data_b[final_data_b$state == "Infected", "counts"]
final_removed_b <- final_data_b[final_data_b$state == "Recovered", "counts"]
final_susceptible_b <- final_data_b[final_data_b$state == "Susceptible", "counts"]

expect_equal(final_infected_b, 0, info = "No agents should remain infected")
expect_equal(final_removed_b, 1, info = "Only the initially infected should recover")
expect_equal(final_susceptible_b, 99, info = "Rest should remain susceptible")

# Test virus with pointer-based parameters -------------------------------------
expect_silent(model_ptr <- ModelSIRCONN(
  name = "Virus Ptr Test",
  n = 50,
  prevalence = 0.1,
  contact_rate = 5,
  transmission_rate = 0.5,
  recovery_rate = 0.3
))

# Add a custom parameter to the model
expect_silent(add_param(model_ptr, "custom_trans", 0.9))

virus_ptr <- virus(
  name = "Pointer Virus",
  prevalence = 0,
  as_proportion = TRUE,
  prob_infecting = 0.1,
  recovery_rate = 0.2
)

# Test pointer-based probability setting
expect_silent(set_prob_infecting_ptr(virus_ptr, model_ptr, "custom_trans"))
expect_silent(set_prob_recovery_ptr(virus_ptr, model_ptr, "Recovery rate"))
expect_silent(add_virus(model_ptr, virus_ptr))
expect_silent(verbose_off(model_ptr))
expect_silent(run(model_ptr, ndays = 5, seed = 789))
