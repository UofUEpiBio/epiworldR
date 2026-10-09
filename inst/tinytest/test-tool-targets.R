# Test just this file: tinytest::run_test_file("inst/tinytest/test-tool-targets.R")

# Two-virus model: virus "A" is created by the model (lineage 0) and virus "B"
# is added afterwards (lineage 1).
make_model <- function() {

  model <- ModelSIRCONN(
    name              = "A",
    n                 = 1000,
    prevalence        = 0.02,
    contact_rate      = 4,
    transmission_rate = 0.5,
    recovery_rate     = 0.2
  )

  verbose_off(model)

  virus_b <- virus(
    name           = "B",
    prevalence     = 20,
    as_proportion  = FALSE,
    prob_infecting = 0.5,
    recovery_rate  = 0.2
  )

  add_virus(model, virus_b)

  list(model = model, virus_b = virus_b)

}

# Everyone gets a tool that fully blocks infection
make_tool <- function() {
  tool(
    name                     = "Perfect vaccine",
    prevalence               = 1,
    as_proportion            = TRUE,
    susceptibility_reduction = 1,
    transmission_reduction   = 0,
    recovery_enhancer        = 0,
    death_reduction          = 0
  )
}

# Number of transmissions (excluding the initial cases) by virus id
n_transmissions <- function(model, virus_id) {
  tn <- get_transmissions(model)
  sum(tn$source >= 0 & tn$virus_id == virus_id)
}

# Lineage ids ------------------------------------------------------------------
m <- make_model()
expect_equal(get_lineage_virus(get_virus(m$model, 0)), 0L)
expect_equal(get_lineage_virus(get_virus(m$model, 1)), 1L)
expect_equal(get_lineage_virus(m$virus_b), 1L) # add_virus() assigns it

virus_c <- virus(
  name = "C", prevalence = 1, as_proportion = FALSE,
  prob_infecting = 0.5, recovery_rate = 0.2
)
expect_true(is.na(get_lineage_virus(virus_c))) # not added to a model yet

# Getters and setters ----------------------------------------------------------
t0 <- make_tool()
expect_identical(get_targets_tool(t0), integer(0)) # default: every virus

expect_silent(add_target_tool(t0, m$virus_b))
expect_identical(get_targets_tool(t0), 1L)

expect_silent(add_target_tool(t0, 0))
expect_identical(get_targets_tool(t0), c(0L, 1L))

expect_silent(set_targets_tool(t0, 5L))
expect_identical(get_targets_tool(t0), 5L)

expect_silent(set_targets_tool(t0, list(get_virus(m$model, 0), m$virus_b)))
expect_identical(get_targets_tool(t0), c(0L, 1L))

expect_silent(set_targets_tool(t0, NULL))
expect_identical(get_targets_tool(t0), integer(0))

expect_silent(set_targets_tool(t0, c(2, 62)))
expect_silent(clear_targets_tool(t0))
expect_identical(get_targets_tool(t0), integer(0))

expect_inherits(add_target_tool(t0, 3), "epiworld_tool")
expect_inherits(set_targets_tool(t0, 3), "epiworld_tool")
expect_inherits(clear_targets_tool(t0), "epiworld_tool")

# Input validation -------------------------------------------------------------
expect_error(add_target_tool(t0, virus_c), "has no lineage id")
expect_error(add_target_tool(t0, 63), "0 to 62")
expect_error(add_target_tool(t0, -1), "0 to 62")
expect_error(add_target_tool(t0, 1.5), "integer lineage ids")
expect_error(add_target_tool(t0, NA), "integer lineage ids")
expect_error(add_target_tool(t0, "A"), "integer lineage ids")
expect_error(add_target_tool(t0, list(m$virus_b, "A")), "epiworld_virus")
expect_error(set_targets_tool(t0, 63), "0 to 62")
expect_error(add_target_tool(m$model, 0), "epiworld_tool")
expect_error(get_targets_tool(m$virus_b), "epiworld_tool")
expect_error(get_lineage_virus(t0), "epiworld_virus")

# A failed set_targets_tool() leaves the targets unchanged
expect_silent(set_targets_tool(t0, 4))
expect_error(set_targets_tool(t0, c(1, 70)))
expect_identical(get_targets_tool(t0), 4L)

# Default: the tool acts on every virus ----------------------------------------
m_all <- make_model()
add_tool(m_all$model, make_tool())
run(m_all$model, ndays = 20, seed = 1231)
expect_equal(n_transmissions(m_all$model, 0), 0)
expect_equal(n_transmissions(m_all$model, 1), 0)

# Targeting A: protects against A but not B ------------------------------------
m_a <- make_model()
t_a <- make_tool()
add_target_tool(t_a, get_virus(m_a$model, 0))
add_tool(m_a$model, t_a)
expect_identical(get_targets_tool(get_tool(m_a$model, 0)), 0L)
run(m_a$model, ndays = 20, seed = 1231)
expect_equal(n_transmissions(m_a$model, 0), 0)
expect_true(n_transmissions(m_a$model, 1) > 0)

# Targeting B (by virus object): protects against B but not A ------------------
m_b <- make_model()
t_b <- make_tool()
add_target_tool(t_b, m_b$virus_b)
add_tool(m_b$model, t_b)
run(m_b$model, ndays = 20, seed = 1231)
expect_true(n_transmissions(m_b$model, 0) > 0)
expect_equal(n_transmissions(m_b$model, 1), 0)

# Targets survive run_multiple() (the model's tools are copied per thread) ----
m_mult <- make_model()
t_mult <- make_tool()
add_target_tool(t_mult, 0)
add_tool(m_mult$model, t_mult)
saver <- make_saver("transmission")
run_multiple(
  m_mult$model, ndays = 20, nsims = 2, seed = 1231, saver = saver,
  nthreads = 1, verbose = FALSE
)
tn_mult <- run_multiple_get_results(m_mult$model, nthreads = 1)$transmission
expect_equal(sum(tn_mult$source >= 0 & tn_mult$virus_id == 0), 0)
expect_true(sum(tn_mult$source >= 0 & tn_mult$virus_id == 1) > 0)
