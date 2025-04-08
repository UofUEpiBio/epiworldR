#' Model Diagram
#'
#' @aliases epiworld_modeldiagram
#' @details
#' Draws model diagrams, such as mermaid diagrams based on transition
#' probability matrices.
#' @returns
#' - The `ModelDiagram` function returns an object of class [epiworld_modeldiagram].
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
#' # Create a ModelDiagram object
#' m_diagram <- ModelDiagram()
#'
#' # Draw mermaid diagrams from model data
#' draw_mermaid_from_data(
#'   m_diagram,
#'   states = get_states(model),
#'   transition_probs = c(get_transition_probability(model))
#' )
#' @export
ModelDiagram <- function() {
  structure(
    ModelDiagram_cpp(),
    class = c("epiworld_modeldiagram")
  )
}

#' @rdname ModelDiagram
#' @param model_diagram ModelDiagram object
#' @param states String vector. List of model states.
#' @param transition_probs Numeric vector. Transition probability matrix
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The ModelDiagram object.
#' @export
draw_mermaid_from_data <- function(
    model_diagram,
    states,
    transition_probs,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_modeldiagram(model_diagram)
  stopifnot_stringvector(states)
  stopifnot_numvector(transition_probs)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  invisible(draw_from_data_cpp(
    model_diagram,
    states,
    transition_probs,
    output_file,
    allow_self_transitions
  ))
}

#' @rdname ModelDiagram
#' @param model_diagram ModelDiagram object
#' @param transitions_file String. Path to file containing the transition probabilities matrix.
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The ModelDiagram object.
#' @export
draw_mermaid_from_file <- function(
    model_diagram,
    transitions_file,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_modeldiagram(model_diagram)
  stopifnot_string(transitions_file)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  invisible(draw_from_file_cpp(
    model_diagram,
    transitions_file,
    output_file,
    allow_self_transitions
  ))
}

#' @rdname ModelDiagram
#' @param model_diagram ModelDiagram object
#' @param transitions_files String vector. List of files containing transition probabilities matrices from multiple model runs.
#' @param output_file String. Path to the output file.
#' @param allow_self_transitions Logical. Whether to allow self-transitions.
#' @returns
#' - The ModelDiagram object.
#' @export
draw_mermaid_from_files <- function(
    model_diagram,
    transitions_files,
    output_file = "",
    allow_self_transitions = FALSE
    ) {
  stopifnot_modeldiagram(model_diagram)
  stopifnot_stringvector(transitions_files)
  stopifnot_string(output_file)
  stopifnot_bool(allow_self_transitions)

  invisible(draw_from_files_cpp(
    model_diagram,
    transitions_files,
    output_file,
    allow_self_transitions
  ))
}
