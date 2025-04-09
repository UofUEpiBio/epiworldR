#' Model Diagram
#'
#' Functions described below are helper functions for drawing
#' diagrams from model data. These generate mermaid diagrams from
#' transition probability matrices which can then be rendered
#' using other packages.
#'
#' @name epiworld-model-diagram
#' @examples
#' # Create and run a model
#' model <- ModelSIRCONN(
#'   name = "A Virus",
#'   n = 10000,
#'   prevalence = .01,
#'   contact_rate = 4.0,
#'   transmission_rate = .5,
#'   recovery_rate = 1.0 / 7.0
#' )
#'
#' verbose_off(model)
#' run(model, ndays = 50, seed = 1912)
#'
#' # Draw mermaid diagrams from model data
#' draw_mermaid_from_data(
#'   states = get_states(model),
#'   transition_probs = c(get_transition_probability(model))
#' )
#' @param states String vector. List of model states.
#' @param transition_probs Numeric vector. Transition probability matrix
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The `draw_mermaid_from_data` function returns the
#' mermaid diagram as a string.
#' @export
draw_mermaid_from_data <- function(
    states,
    transition_probs,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_stringvector(states)
  stopifnot_numvector(transition_probs)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  diagram <- capture.output(draw_from_data_cpp(
    states,
    transition_probs,
    output_file,
    allow_self_transitions
  ))

  return(paste(diagram, collapse = "\n"))
}

#' @rdname epiworld-model-diagram
#' @param transitions_file String. Path to file containing the transition probabilities matrix.
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The `draw_mermaid_from_file` function returns the
#' mermaid diagram as a string.
#' @export
draw_mermaid_from_file <- function(
    transitions_file,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_string(transitions_file)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  diagram <- capture.output(draw_from_file_cpp(
    transitions_file,
    output_file,
    allow_self_transitions
  ))

  return(paste(diagram, collapse = "\n"))
}

#' @rdname epiworld-model-diagram
#' @param transitions_files String vector. List of files containing transition probabilities matrices from multiple model runs.
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The `draw_mermaid_from_files` function returns the
#' mermaid diagram as a string.
#' @export
draw_mermaid_from_files <- function(
    transitions_files,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_stringvector(transitions_files)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  diagram <- capture.output(draw_from_files_cpp(
    transitions_files,
    output_file,
    allow_self_transitions
  ))

  return(paste(diagram, collapse = "\n"))
}
