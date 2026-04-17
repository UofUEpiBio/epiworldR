# Test just this file: tinytest::run_test_file("inst/tinytest/test-model-builder.R")

# Build a SEIR model using the model_builder framework and compare it with the
# built-in ModelSEIR implementation.

n                 <- 1000L
prevalence        <- 0.01
trans_rate        <- 0.5
incub_days        <- 7.0
recov_rate        <- 0.1
ndays             <- 50L
seed              <- 1912L
# Maximum allowed absolute deviation (as proportion of n) between the two
# implementations.  The two models use slightly different RNG call patterns
# (roulette vs. direct Bernoulli comparison) so exact equality is not expected.
comparison_tol <- 0.15

# ---- Reference: built-in ModelSEIR on a small-world network -----------------

seir_ref <- ModelSEIR(
  name              = "COVID-19",
  prevalence        = prevalence,
  transmission_rate = trans_rate,
  incubation_days   = incub_days,
  recovery_rate     = recov_rate
)
agents_smallworld(seir_ref, n = n, k = 5L, d = FALSE, p = 0.01)
verbose_off(seir_ref)
run(seir_ref, ndays = ndays, seed = seed)

# ---- model_builder SEIR (S-E-I-R with identical compartments) ---------------

# Create a blank model and add parameters matching ModelSEIR
seir_mb <- new_model()
add_param(seir_mb, "Transmission rate", trans_rate)
add_param(seir_mb, "Incubation rate",   1.0 / incub_days)
add_param(seir_mb, "Recovery rate",     recov_rate)

# Update functions:
#   S (0): susceptible sampler; exclude Removed (3) so that agents who have
#          already recovered (but still carry the virus pointer) do not infect.
#   E (1): transitions to Infected (2) at the incubation rate
#   I (2): transitions to Removed (3) at the recovery rate
#   R (3): no update (absorbing)
update_s <- update_fun_susceptible(exclude = 3L)
update_e <- update_fun_rate(param_names = "Incubation rate", target_states = 2L)
update_i <- update_fun_rate(param_names = "Recovery rate",   target_states = 3L)

add_state(seir_mb, "Susceptible", update_s)
add_state(seir_mb, "Exposed",     update_e)
add_state(seir_mb, "Infected",    update_i)
add_state(seir_mb, "Removed",     NULL)

# Virus: agents start in Exposed (1) when first infected; end/removed → Removed (3)
covid <- virus(
  name          = "COVID-19",
  prevalence    = prevalence,
  as_proportion = TRUE,
  prob_infecting = trans_rate,
  recovery_rate  = recov_rate
)
virus_set_state(covid, init = 1L, end = 3L, removed = 3L)
set_prob_infecting_ptr(covid, seir_mb, "Transmission rate")
set_distribution_virus(covid, distribute_virus_randomly(prevalence, TRUE))
add_virus(seir_mb, covid)

agents_smallworld(seir_mb, n = n, k = 5L, d = FALSE, p = 0.01)
verbose_off(seir_mb)

expect_silent(run(seir_mb, ndays = ndays, seed = seed))

# ---- Structural checks -------------------------------------------------------

# Both models should have the same four state names in the same order
expect_equal(get_states(seir_mb), get_states(seir_ref))

# Population must be conserved at every time step
hist_mb <- get_hist_total(seir_mb)
totals_per_day <- tapply(hist_mb$counts, hist_mb$date, sum)
expect_true(all(totals_per_day == n))

# ---- Epidemiological comparison ---------------------------------------------

# Extract final-day counts for each model
final_ref <- get_hist_total(seir_ref)
final_ref <- final_ref[final_ref$date == ndays, c("state", "counts")]

final_mb <- hist_mb[hist_mb$date == ndays, c("state", "counts")]

# The total removed (immune) should be within 15 % of the reference model.
# (The two implementations use slightly different RNG call patterns, so exact
# equality is not expected, but both should capture the same epidemic dynamics.)
removed_ref <- final_ref$counts[final_ref$state == "Removed"]
removed_mb  <- final_mb$counts[final_mb$state  == "Removed"]

expect_true(
  abs(removed_mb - removed_ref) / n < comparison_tol,
  info = sprintf(
    "Final Removed counts differ by more than 15%%: ref=%d, mb=%d (n=%d)",
    removed_ref, removed_mb, n
  )
)
