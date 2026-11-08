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
# The registry: models by name -------------------------------------------------
models <- list_models()
expect_true(all(
  c("SIS", "SIR", "SEIR", "SISD", "SIRD", "SEIRD", "SEIRQuarantine", "SURV") %in%
    models$id
))

info <- get_model_info("SEIR")
expect_identical(
  info$params$arg,
  c("n", "prevalence", "transmission_prob", "latent_period", "recovery_rate")
)
expect_equal(info$params$name[3], "Transmission probability")

# A named constructor is make_model() with the same arguments: the same run
by_name <- make_model(
  "SEIR",
  list(
    n = 2000, prevalence = .01, transmission_prob = .3,
    latent_period = 3, recovery_rate = .1
  ),
  virus_name = "flu"
)
shim <- ModelSEIR("flu", 2000, .01, .3, 3, .1)
expect_identical(class(by_name), class(shim))
for (m in list(by_name, shim)) {
  add_layer(m, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
  verbose_off(m)
  run(m, ndays = 30, seed = 5)
}
expect_identical(get_hist_total(by_name), get_hist_total(shim))

# Arguments left out take their defaults; values outside the limits fail
defaults <- make_model("SIR")
expect_equal(size(defaults), get_model_info("SIR")$params$value[1])
expect_error(make_model("SIR", list(transmission_prob = 1.5)), "transmission_prob")
expect_error(make_model("SIR", list(bad_arg = 1)), "bad_arg")
expect_error(make_model("XYZ"), "SIR")
expect_error(make_model("SIR", list(.1)), "named")
