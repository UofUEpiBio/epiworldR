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
model_a <- ModelSIR(
  name = "Virus Test (a)",
  n = 100,
  prevalence = 0.01,
  transmission_prob = 1.0,
  recovery_rate = 0.0
)
add_param(model_a, "Contact rate", 10)
add_layer(
  model_a,
  layer_mixing("everyone", groups = rep(0, 100), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

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
model_b <- ModelSIR(
  name = "Virus Test (b)",
  n = 100,
  prevalence = 0.01,
  transmission_prob = 0.0,
  recovery_rate = 1.0
)
add_param(model_b, "Contact rate", 10)
add_layer(
  model_b,
  layer_mixing("everyone", groups = rep(0, 100), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

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
model_ptr <- ModelSIR(
  name = "Virus Ptr Test",
  n = 50,
  prevalence = 0.1,
  transmission_prob = 0.5,
  recovery_rate = 0.3
)
add_param(model_ptr, "Contact rate", 5)
add_layer(
  model_ptr,
  layer_mixing("everyone", groups = rep(0, 50), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

# Add custom parameters to the model
expect_silent(add_param(model_ptr, "custom_trans", 0.9))
expect_silent(add_param(model_ptr, "custom_death", 0.05))
expect_silent(add_param(model_ptr, "custom_incubation", 3.0))

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
expect_silent(set_prob_death_ptr(virus_ptr, model_ptr, "custom_death"))
expect_silent(set_incubation_ptr(virus_ptr, model_ptr, "custom_incubation"))

# Verify parameters can be retrieved
expect_equal(get_param(model_ptr, "custom_trans"), 0.9)
expect_equal(get_param(model_ptr, "custom_death"), 0.05)
expect_equal(get_param(model_ptr, "custom_incubation"), 3.0)

expect_silent(add_virus(model_ptr, virus_ptr))
expect_silent(verbose_off(model_ptr))
expect_silent(run(model_ptr, ndays = 5, seed = 789))

# Adding a second virus to models on a mixing layer ----------------------------
cmat <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2) * 5
groups <- rep(0:1, each = 250)
mixing_models <- list(
  ModelSIR("A virus", 500, 0.01, 0.2, 1 / 7),
  ModelSEIR("A virus", 500, 0.01, 0.2, 3, 1 / 7)
)
for (m in mixing_models) {
  add_layer(m, layer_mixing("mixing", groups, contact_matrix = cmat))
  other <- virus("Another virus", 0.02, TRUE, 0.2, recovery_rate = 1 / 7)
  expect_silent(add_virus(m, other))
  verbose_off(m)
  expect_silent(run(m, ndays = 10, seed = 1))
}

quarantine_model <- ModelSEIRQuarantine(
  name = "A virus", n = 500, prevalence = 0.01, transmission_prob = 0.2,
  recovery_rate = 1 / 7, latent_period = 3,
  hospitalization_rate = 0.05, hospitalization_period = 7,
  detection_delay = 3, quarantine_period = 14, quarantine_willingness = 0.8,
  isolation_willingness = 0.5, isolation_period = 7,
  contact_tracing_prob = 0.7, contact_tracing_window = 3
)
add_layer(quarantine_model, layer_mixing("mixing", groups, contact_matrix = cmat))
expect_silent(add_virus(
  quarantine_model,
  virus("Another virus", 0.02, TRUE, 0.2, recovery_rate = 1 / 7)
))
