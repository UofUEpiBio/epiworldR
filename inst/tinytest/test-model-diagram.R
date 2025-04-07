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

# Check draw_mermaid_from_data
expect_stdout(draw_mermaid_from_data(
    m_diagram,
    states = get_states(model),
    transition_probs = c(get_transition_probability(model))
))
expect_stdout(draw_mermaid_from_data(
  m_diagram,
  states = get_states(model),
  transition_probs = c(get_transition_probability(model)),
  allow_self_transitions = TRUE
))
expect_silent(draw_mermaid_from_data(
  m_diagram,
  states = get_states(model),
  transition_probs = c(get_transition_probability(model)),
  output_file = "mermaid_diagram.txt"
))
expect_true(file.exists("mermaid_diagram.txt"))
file.remove("mermaid_diagram.txt")
expect_false(file.exists("mermaid_diagram.txt"))

# Check draw_mermaid_from_file
test_file <- "test-model-diagram-files/single-file-transitions.txt"
expect_stdout(draw_mermaid_from_file(
    m_diagram,
    transitions_file = test_file
))
expect_stdout(draw_mermaid_from_file(
    m_diagram,
    transitions_file = test_file,
    allow_self_transitions = TRUE
))
expect_silent(draw_mermaid_from_file(
    m_diagram,
    transitions_file = test_file,
    output_file = "mermaid_diagram.txt"
))
expect_true(file.exists("mermaid_diagram.txt"))
file.remove("mermaid_diagram.txt")
expect_false(file.exists("mermaid_diagram.txt"))

# Check draw_mermaid_from_files
test_files <- paste0("test-model-diagram-files/multiple-files/", list.files("test-model-diagram-files/multiple-files/"))
expect_stdout(draw_mermaid_from_files(
    m_diagram,
    transitions_files = test_files
))
expect_stdout(draw_mermaid_from_files(
    m_diagram,
    transitions_files = test_files,
    allow_self_transitions = TRUE
))
expect_silent(draw_mermaid_from_files(
    m_diagram,
    transitions_files = test_files,
    output_file = "mermaid_diagram.txt"
))
expect_true(file.exists("mermaid_diagram.txt"))
file.remove("mermaid_diagram.txt")
expect_false(file.exists("mermaid_diagram.txt"))

# Check functions fail with invalid inputs -------------------------------------
# Check draw_mermaid_from_data
expect_error(draw_mermaid_from_data(
    NULL,
    states = get_states(model),
    transition_probs = c(get_transition_probability(model))
), "must be an object of class 'epiworld_modeldiagram'")
expect_error(draw_mermaid_from_data(
    m_diagram,
    states = NA,
    transition_probs = c(get_transition_probability(model))
), "must be a string vector")
expect_error(draw_mermaid_from_data(
    m_diagram,
    states = get_states(model),
    transition_probs = NA
), "must be a numeric vector")
expect_error(draw_mermaid_from_data(
    m_diagram,
    states = get_states(model),
    transition_probs = c(get_transition_probability(model)),
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_data(
    m_diagram,
    states = get_states(model),
    transition_probs = c(get_transition_probability(model)),
    allow_self_transitions = "Not a Bool"
), "must be a boolean")

# Check draw_mermaid_from_file
expect_error(draw_mermaid_from_file(
    NA,
    transitions_file = test_file
), "must be an object of class 'epiworld_modeldiagram'")
expect_error(draw_mermaid_from_file(
    m_diagram,
    transitions_file = NA
), "must be a string")
expect_error(draw_mermaid_from_file(
    m_diagram,
    transitions_file = test_file,
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_file(
    m_diagram,
    transitions_file = test_file,
    allow_self_transitions = "Not a Bool"
), "must be a boolean")

# Check draw_mermaid_from_files
expect_error(draw_mermaid_from_files(
    NA,
    transitions_files = test_files
), "must be an object of class 'epiworld_modeldiagram'")
expect_error(draw_mermaid_from_files(
    m_diagram,
    transitions_files = NA
), "must be a string vector")
expect_error(draw_mermaid_from_files(
    m_diagram,
    transitions_files = test_files,
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_files(
    m_diagram,
    transitions_files = test_files,
    allow_self_transitions = "Not a Bool"
), "must be a boolean")