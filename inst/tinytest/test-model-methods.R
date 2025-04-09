# Test just this file: tinytest::run_test_file("inst/tinytest/test-model-methods.R")

# Setup for tests --------------------------------------------------------------
model <- ModelSIRCONN(
  name = "A Virus",
  n = 10000,
  prevalence = .01,
  contact_rate = 4.0,
  transmission_rate = .5,
  recovery_rate = 1.0/7.0
)

verbose_off(model)
run(model, ndays = 50, seed = 1912)

# Check draw_mermaid succeeds with valid inputs --------------------------------
expect_silent(draw_mermaid(model))
expect_silent(draw_mermaid(model, allow_self_transitions = TRUE))
expect_silent(draw_mermaid(model, output_file = "mermaid_diagram.txt"))
expect_true(file.exists("mermaid_diagram.txt"))
file.remove("mermaid_diagram.txt")
expect_false(file.exists("mermaid_diagram.txt"))

# Check functions fail with invalid inputs -------------------------------------

# Check draw_mermaid function
expect_error(draw_mermaid(NULL), "must be of class 'epiworld_model'.")
expect_error(draw_mermaid(model, output_file = NULL), "must be a string")
expect_error(draw_mermaid(model, allow_self_transitions = NA), "must not be NA")
expect_error(draw_mermaid(model, allow_self_transitions = "NOT_BOOL"), "must be a boolean")