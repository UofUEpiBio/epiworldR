# Test just this file: tinytest::run_test_file("inst/tinytest/test-clone-epi.R")

make_model <- function() {
  model <- ModelSIRCONN(
    name              = "A",
    n                 = 1000,
    prevalence        = 0.05,
    contact_rate      = 4,
    transmission_rate = 0.5,
    recovery_rate     = 0.2
  )
  verbose_off(model)
  model
}

perfect_vaccine <- function() {
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

final_recovered <- function(model) {
  hist <- get_hist_total(model)
  hist[hist$date == max(hist$date) & hist$state == "Recovered", "counts"]
}

# ------------------------------------------------------------------------------
# Tools
# ------------------------------------------------------------------------------
vax  <- perfect_vaccine()
vax2 <- clone_epi(vax)

expect_inherits(vax2, "epiworld_tool")
expect_equal(get_name_tool(vax2), "Perfect vaccine")

# Changing the copy does not change the original (unlike `vax2 <- vax`)
set_name_tool(vax2, "Useless vaccine")
set_susceptibility_reduction(vax2, 0)
expect_equal(get_name_tool(vax), "Perfect vaccine")
expect_equal(get_name_tool(vax2), "Useless vaccine")

# Everyone gets the original (perfect) vaccine: only the initial cases
model <- make_model()
add_tool(model, vax)
run(model, ndays = 50, seed = 123)
expect_equal(final_recovered(model), 50)

# Everyone gets the copy (useless) vaccine: the outbreak spreads
model <- make_model()
add_tool(model, vax2)
run(model, ndays = 50, seed = 123)
expect_true(final_recovered(model) > 100)

# The model records the original and the copy as different tools (fresh
# objects, since adding a tool to a model stamps it with that model's id)
vax  <- perfect_vaccine()
vax2 <- clone_epi(vax)
set_name_tool(vax2, "Useless vaccine")
model <- make_model()
add_tool(model, vax)
add_tool(model, vax2)
expect_equal(get_n_tools(model), 2L)
expect_equal(get_name_tool(get_tool(model, 0)), "Perfect vaccine")
expect_equal(get_name_tool(get_tool(model, 1)), "Useless vaccine")
run(model, ndays = 10, seed = 123)
expect_equal(sort(unique(get_hist_tool(model)$tool_id)), c(0L, 1L))

# Cloning a tool that is already in a model also gives a new tool
model <- make_model()
add_tool(model, vax)
vax3 <- clone_epi(get_tool(model, 0))
set_name_tool(vax3, "Copy")
add_tool(model, vax3)
run(model, ndays = 10, seed = 123)
expect_equal(get_name_tool(get_tool(model, 0)), "Perfect vaccine")
expect_equal(sort(unique(get_hist_tool(model)$tool_id)), c(0L, 1L))

# ------------------------------------------------------------------------------
# Viruses
# ------------------------------------------------------------------------------
model <- make_model()
virus_b <- clone_epi(get_virus(model, 0))
expect_inherits(virus_b, "epiworld_virus")
set_name_virus(virus_b, "B")
add_virus(model, virus_b)

expect_equal(get_n_viruses(model), 2L)
expect_equal(get_name_virus(get_virus(model, 0)), "A")
expect_equal(get_name_virus(get_virus(model, 1)), "B")

# The copy founds its own lineage
run(model, ndays = 10, seed = 123)
expect_equal(get_lineage_virus(get_virus(model, 0)), 0L)
expect_equal(get_lineage_virus(get_virus(model, 1)), 1L)

# ------------------------------------------------------------------------------
# Models and other objects
# ------------------------------------------------------------------------------
model  <- make_model()
model2 <- clone_epi(model)
expect_inherits(model2, class(model))

# The copy is independent of the original
name0 <- get_name(model)
set_name(model2, "Copy")
expect_equal(get_name(model2), "Copy")
expect_equal(get_name(model), name0)

# clone_model() is deprecated in favor of clone_epi()
expect_warning(model3 <- clone_model(model), "clone_epi")
expect_inherits(model3, class(model))
expect_error(clone_epi(1), "no method")
