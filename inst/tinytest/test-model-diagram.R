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

# Check draw functions succeed with valid inputs -------------------------------

test_output_file <- "mermaid_diagram.txt"

# Check draw_mermaid_from_data
expect_false(identical(
    draw_mermaid_from_data(
        states = get_states(model),
        transition_probs = c(get_transition_probability(model))
    ),
    draw_mermaid_from_data(
        states = get_states(model),
        transition_probs = c(get_transition_probability(model)),
        allow_self_transitions = TRUE
    )
))

expect_identical(
    draw_mermaid_from_data(
        states = get_states(model),
        transition_probs = c(get_transition_probability(model))
    ),
    draw_mermaid_from_data(
        states = get_states(model),
        transition_probs = c(get_transition_probability(model)),
        output_file = test_output_file
    )
)

expect_true(file.exists(test_output_file))
file.remove(test_output_file)
expect_false(file.exists(test_output_file))

# Check draw_mermaid_from_matrix
expect_false(identical(
    draw_mermaid_from_matrix(
        transition_matrix = get_transition_probability(model)
    ),
    draw_mermaid_from_matrix(
        transition_matrix = get_transition_probability(model),
        allow_self_transitions = TRUE
    )    
))

expect_identical(
    draw_mermaid_from_matrix(
        transition_matrix = get_transition_probability(model)
    ),
    draw_mermaid_from_matrix(
        transition_matrix = get_transition_probability(model),
        output_file = test_output_file
    )
)

expect_true(file.exists(test_output_file))
file.remove(test_output_file)
expect_false(file.exists(test_output_file))

# Check draw_mermaid_from_file
test_file <- "test-model-diagram-files/single-file-transitions.txt"
expect_false(identical(
    draw_mermaid_from_file(
        transitions_file = test_file,
    ),
    draw_mermaid_from_file(
        transitions_file = test_file,
        allow_self_transitions = TRUE
    )
))

expect_identical(
    draw_mermaid_from_file(
        transitions_file = test_file
    ), 
    draw_mermaid_from_file(
        transitions_file = test_file,
        output_file = test_output_file
    )
)
expect_true(file.exists(test_output_file))
file.remove(test_output_file)
expect_false(file.exists(test_output_file))

# Check draw_mermaid_from_files
test_files <- paste0("test-model-diagram-files/multiple-files/", list.files("test-model-diagram-files/multiple-files/"))
expect_false(identical(
    draw_mermaid_from_files(
        transitions_files = test_files
    ),
    draw_mermaid_from_files(
        transitions_files = test_files,
        allow_self_transitions = TRUE
    )
))

expect_identical(
    draw_mermaid_from_files(
        transitions_files = test_files
    ),
    draw_mermaid_from_files(
        transitions_files = test_files,
        output_file = test_output_file
    )
)

expect_true(file.exists(test_output_file))
file.remove(test_output_file)
expect_false(file.exists(test_output_file))

# Check functions fail with invalid inputs -------------------------------------
# Check draw_mermaid_from_data
expect_error(draw_mermaid_from_data(
    states = NA,
    transition_probs = c(get_transition_probability(model))
), "must be a string vector")
expect_error(draw_mermaid_from_data(
    states = get_states(model),
    transition_probs = NA
), "must be a numeric vector")
expect_error(draw_mermaid_from_data(
    states = get_states(model),
    transition_probs = c(get_transition_probability(model)),
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_data(
    states = get_states(model),
    transition_probs = c(get_transition_probability(model)),
    allow_self_transitions = "Not a Bool"
), "must be a boolean")

# Check draw_mermaid_from_matrix
expect_error(draw_mermaid_from_matrix(
    transition_matrix = NA,
), "must not contain NA values")

mismatch_rowcol_names_matrix <- matrix(1:9, nrow = 3, dimnames = list(c("X","Y","Z"), c("A","B","C")))
expect_error(draw_mermaid_from_matrix(
  transition_matrix = mismatch_rowcol_names_matrix ,
), "must have the same row and column names")

nonsquare_matrix <- matrix(1:8, nrow = 2, ncol = 4)
expect_error(draw_mermaid_from_matrix(
  transition_matrix = nonsquare_matrix,
), "matrix must be square")

expect_error(draw_mermaid_from_matrix(
    transition_matrix = get_transition_probability(model),
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_matrix(
    transition_matrix = get_transition_probability(model),
    allow_self_transitions = "Not a Bool"
), "must be a boolean")

# Check draw_mermaid_from_file
expect_error(draw_mermaid_from_file(
    transitions_file = NA
), "must be a string")
expect_error(draw_mermaid_from_file(
    transitions_file = test_file,
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_file(
    transitions_file = test_file,
    allow_self_transitions = "Not a Bool"
), "must be a boolean")

# Check draw_mermaid_from_files
expect_error(draw_mermaid_from_files(
    transitions_files = NA
), "must be a string vector")
expect_error(draw_mermaid_from_files(
    transitions_files = test_files,
    output_file = NA
), "must be a string")
expect_error(draw_mermaid_from_files(
    transitions_files = test_files,
    allow_self_transitions = "Not a Bool"
), "must be a boolean")