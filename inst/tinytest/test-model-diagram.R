# Test just this file: tinytest::run_test_file("inst/tinytest/test-model-diagram.R")

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

# Check initialization ---------------------------------------------------------
expect_silent(m_diagram <- ModelDiagram())
expect_inherits(m_diagram, "epiworld_modeldiagram")
expect_length(class(m_diagram), 1)

# Check draw functions succeed with valid inputs -------------------------------

# Check draw_from_data


# Check draw_from_file

# Check draw_from_files

# Check functions fail with invalid inputs -------------------------------------