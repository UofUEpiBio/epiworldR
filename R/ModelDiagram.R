#' Model Diagram
#'
#' Functions described here are helper functions for drawing
#' diagrams from model data. These generate mermaid diagrams from
#' transition probability matrices which can then be rendered
#' using other packages.
#'
#' @name epiworld-model-diagram
#' @concept model-utility-functions
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
#' @param output_file String. Optional path to a file. If provided, the diagram will be written to the file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions, defaults to FALSE.
#' @returns
#' - The `draw_mermaid_from_data` function returns the
#' mermaid diagram as a string.
#' @export
#' @importFrom utils capture.output
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

  if (output_file != "") {
    draw_from_data_cpp(
      states,
      transition_probs,
      output_file,
      allow_self_transitions
    )

    message("Diagram written to ", output_file)

    diagram <- readChar(output_file, file.info(output_file)$size)
    return(diagram)
  } else {
    diagram <- capture.output(draw_from_data_cpp(
      states,
      transition_probs,
      output_file,
      allow_self_transitions
    ))

    return(paste(diagram, collapse = "\n"))
  }
}

#' @rdname epiworld-model-diagram
#' @param transition_matrix Square matrix. Contains states and transition probabilities.
#' @returns
#' - The `draw_mermaid_from_matrix` function returns the
#' mermaid diagram as a string.
#' @export
draw_mermaid_from_matrix <- function(
  transition_matrix,
  output_file = "",
  allow_self_transitions = FALSE
) {

  stopifany_na(transition_matrix)

  if (nrow(transition_matrix) != ncol(transition_matrix)) {
    stop(paste0(
      "Transition matrix must be square, but instead has dimensions: [",
      nrow(transition_matrix), "x", ncol(transition_matrix), "]"))
  }

  if (!identical(colnames(transition_matrix), rownames(transition_matrix))) {
    stop(paste0(
      "Transition matrix must have the same row and column names, but instead has row names: [",
      paste(rownames(transition_matrix), collapse = ", "), "] and column names: [",
      paste(colnames(transition_matrix), collapse = ", "), "]"))
  }

  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  m_states <- colnames(transition_matrix)
  t_probs <- c(transition_matrix)

  draw_mermaid_from_data(
    states = m_states,
    transition_probs = t_probs,
    output_file = output_file,
    allow_self_transitions = allow_self_transitions
  )
}

#' @rdname epiworld-model-diagram
#' @param transitions_file String. Path to file containing the transition probabilities matrix.
#' @returns
#' - The `draw_mermaid_from_file` function returns the
#' mermaid diagram as a string.
#' @export
#' @importFrom utils capture.output
draw_mermaid_from_file <- function(
  transitions_file,
  output_file = "",
  allow_self_transitions = FALSE
) {
  stopifnot_string(transitions_file)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  if (output_file != "") {
    draw_from_file_cpp(
      transitions_file,
      output_file,
      allow_self_transitions
    )

    message("Diagram written to ", output_file)

    diagram <- readChar(output_file, file.info(output_file)$size)
    return(diagram)
  } else {
    diagram <- capture.output(draw_from_file_cpp(
      transitions_file,
      output_file,
      allow_self_transitions
    ))

    return(paste(diagram, collapse = "\n"))
  }
}

#' @rdname epiworld-model-diagram
#' @param transitions_files String vector. List of files containing transition probabilities matrices from multiple model runs.
#' @returns
#' - The `draw_mermaid_from_files` function returns the
#' mermaid diagram as a string.
#' @export
#' @importFrom utils capture.output
draw_mermaid_from_files <- function(
  transitions_files,
  output_file = "",
  allow_self_transitions = FALSE
) {
  stopifnot_stringvector(transitions_files)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  if (output_file != "") {
    draw_from_files_cpp(
      transitions_files,
      output_file,
      allow_self_transitions
    )

    message("Diagram written to ", output_file)

    diagram <- readChar(output_file, file.info(output_file)$size)
    return(diagram)
  } else {
    diagram <- capture.output(draw_from_files_cpp(
      transitions_files,
      output_file,
      allow_self_transitions
    ))

    return(paste(diagram, collapse = "\n"))
  }
}
