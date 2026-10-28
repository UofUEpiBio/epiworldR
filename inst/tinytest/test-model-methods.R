# Test just this file: tinytest::run_test_file("inst/tinytest/test-model-methods.R")

# Setup for tests --------------------------------------------------------------
model <- ModelSIR(
  name = "A Virus",
  n = 10000,
  prevalence = .01,
  transmission_prob = .5,
  recovery_rate = 1.0/7.0
)
add_param(model, "Contact rate", 4.0)
add_layer(
  model,
  layer_mixing("everyone", groups = rep(0, 10000), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

verbose_off(model)
run(model, ndays = 50, seed = 1912)

# Check draw_mermaid succeeds with valid inputs --------------------------------
expect_silent(md_basic <- draw_mermaid(model))
expect_silent(md_self <- draw_mermaid(model, allow_self_transitions = TRUE))
expect_message(md_with_output <- draw_mermaid(
  model,
  output_file = tempfile(fileext = ".md")
), "Diagram written")

# Check functions fail with invalid inputs -------------------------------------

# Check draw_mermaid function
expect_error(draw_mermaid(NULL), "must be of class 'epiworld_model'.")
expect_error(draw_mermaid(model, output_file = NULL), "must be a string")
expect_error(draw_mermaid(model, allow_self_transitions = NA), "must not be NA")
expect_error(draw_mermaid(model, allow_self_transitions = "NOT_BOOL"), "must be a boolean")