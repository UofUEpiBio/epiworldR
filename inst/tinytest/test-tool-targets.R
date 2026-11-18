# Test just this file: tinytest::run_test_file("inst/tinytest/test-tool-targets.R")

# Two diseases: "A" is created by the model; "B", built from scratch, is a
# second disease. An agent carries one infection at a time.
two_disease_model <- function() {

  model <- ModelSIR(
    name              = "A",
    n                 = 1000,
    prevalence        = 0.02,
    transmission_prob = 0.5,
    recovery_rate     = 0.2
  )
  add_param(model, "Contact rate", 4)
  add_layer(
    model,
    layer_mixing("everyone", groups = rep(0, 1000), contact_matrix = 1) |>
      set_contact_rate("Contact rate")
  )
  verbose_off(model)

  b <- virus("B") |>
    add_stage("Susceptible") |>
    add_stage("Infected", infectious = TRUE) |>
    add_stage("Recovered") |>
    add_infection("Susceptible", "Infected", rate = 0.5) |>
    add_rate("Infected", "Recovered", rate = 0.2)
  add_virus(model, b, assign_count(20))

  model

}

# A tool that fully blocks infection
make_tool <- function() {
  tool(
    name                     = "Perfect vaccine",
    susceptibility_reduction = 1,
    transmission_reduction   = 0,
    recovery_enhancer        = 0,
    death_reduction          = 0
  )
}

# Number of transmissions (excluding the initial cases) of a disease
n_transmissions <- function(model, disease) {
  tn <- get_transmissions(model, disease = disease)
  sum(tn$source >= 0)
}

# Getters and setters ----------------------------------------------------------
t0 <- make_tool()
expect_identical(get_targets_tool(t0), character(0)) # default: every disease

expect_silent(add_target_tool(t0, "A"))
expect_silent(add_target_tool(t0, "B"))
expect_identical(get_targets_tool(t0), c("A", "B"))

expect_silent(set_targets_tool(t0, "B"))
expect_identical(get_targets_tool(t0), "B")

expect_silent(set_targets_tool(t0, NULL))
expect_identical(get_targets_tool(t0), character(0))

expect_silent(set_targets_tool(t0, c("A", "B")))
expect_silent(clear_targets_tool(t0))
expect_identical(get_targets_tool(t0), character(0))

expect_inherits(add_target_tool(t0, "A"), "epiworld_tool")
expect_inherits(set_targets_tool(t0, "A"), "epiworld_tool")
expect_inherits(clear_targets_tool(t0), "epiworld_tool")

# Input validation: diseases are names, not viruses or ids
m <- two_disease_model()
expect_error(add_target_tool(t0, 0), "disease's name")
expect_error(add_target_tool(t0, get_virus(m, 0)), "disease's name")
expect_error(add_target_tool(t0, c("A", "B")), "disease's name")
expect_error(set_targets_tool(t0, 1:2), "disease names")
expect_error(add_target_tool(m, "A"), "epiworld_tool")
expect_error(get_targets_tool(get_virus(m, 0)), "epiworld_tool")

# The defunct lineage functions say what to use
expect_error(get_lineage_virus(get_virus(m, 0)), "add_target_tool")
expect_error(add_virus_target(get_layer(m, 0), get_virus(m, 0)), "add_target_layer")

# With several diseases, a tool must name its diseases -------------------------
m_none <- two_disease_model()
add_tool(m_none, make_tool(), assign_share(1))
expect_error(run(m_none, ndays = 5, seed = 1), "target")

# A disease the model does not have is an error when the model runs
m_bad <- two_disease_model()
t_bad <- make_tool()
add_target_tool(t_bad, "C")
add_tool(m_bad, t_bad, assign_share(1))
expect_error(run(m_bad, ndays = 5, seed = 1), "C")

# Targeting A: protects against A but not B ------------------------------------
m_a <- two_disease_model()
t_a <- make_tool()
add_target_tool(t_a, "A")
add_tool(m_a, t_a, assign_share(1))
expect_identical(get_targets_tool(get_tool(m_a, 0)), "A")
run(m_a, ndays = 20, seed = 1231)
expect_equal(n_transmissions(m_a, "A"), 0)
expect_true(n_transmissions(m_a, "B") > 0)

# Its history is kept in A only
expect_equal(nrow(get_hist_tool(m_a, disease = "B")), 0L)
expect_true(nrow(get_hist_tool(m_a, disease = "A")) > 0L)

# Targeting B: protects against B but not A ------------------------------------
m_b <- two_disease_model()
t_b <- make_tool()
add_target_tool(t_b, "B")
add_tool(m_b, t_b, assign_share(1))
run(m_b, ndays = 20, seed = 1231)
expect_true(n_transmissions(m_b, "A") > 0)
expect_equal(n_transmissions(m_b, "B"), 0)

# With one disease, a tool without targets acts on it --------------------------
m_one <- ModelSIR(
  name              = "A",
  n                 = 1000,
  prevalence        = 0.02,
  transmission_prob = 0.5,
  recovery_rate     = 0.2
)
add_layer(m_one, layer_network("contacts", graph_smallworld(k = 6, p = .1)))
verbose_off(m_one)
add_tool(m_one, make_tool(), assign_share(1))
run(m_one, ndays = 20, seed = 1231)
expect_equal(n_transmissions(m_one, "A"), 0)

# Targets survive run_multiple() (the model's tools are copied per thread) ----
m_mult <- two_disease_model()
t_mult <- make_tool()
add_target_tool(t_mult, "A")
add_tool(m_mult, t_mult, assign_share(1))
saver <- make_saver("transmission")
run_multiple(
  m_mult, ndays = 20, nsims = 2, seed = 1231, saver = saver,
  nthreads = 1, verbose = FALSE
)
tn_mult <- run_multiple_get_results(m_mult, nthreads = 1)$transmission
expect_equal(sum(tn_mult$source >= 0 & tn_mult$disease == "A"), 0)
expect_true(sum(tn_mult$source >= 0 & tn_mult$disease == "B") > 0)
