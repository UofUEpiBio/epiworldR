#' Tools in epiworld
#' @param m Model
#' @param name Name of the tool
#' @param susceptibility_reduction Numeric. Proportion it reduces susceptibility.
#' @param transmission_reduction Numeric. Proportion it reduces transmission.
#' @param recovery_enhancer Numeric. Proportion it improves recovery.
#' @param death_reduction Numeric. Proportion it reduces probability of death.e
#' @param tool_pos Positive integer. Index of the tool's position in the model.
#' @examples 
#' epitool <- tool(
#'   name = "Vaccine",
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5, 
#'   death_reduction = .9
#' )
#' @export
tool <- function(
    name,
    susceptibility_reduction,
    transmission_reduction,
    recovery_enhancer,
    death_reduction
) {

  structure(
    tool_cpp(
      name,
      susceptibility_reduction,
      transmission_reduction,
      recovery_enhancer,
      death_reduction
    ),
    class = "epiworld_tool"
  )
    
}

#' @export
#' @param t An object of class `epiworld_tool`
#' @param prevalence In the case of `add_tool`, a proportion, otherwise, an integer.
#' @rdname tool
add_tool <- function(m, t, prevalence) UseMethod("add_tool")

#' @export
add_tool.epiworld_model <- function(m, t, prevalence) {
  add_tool_cpp(m, t, prevalence)
  invisible(m)
}

#' @export
#' @rdname tool
add_tool_n <- function(m, t, prevalence) UseMethod("add_tool_n")

#' @export
add_tool_n.epiworld_model <- function(m, t, prevalence) {
  add_tool_n_cpp(m, t, prevalence)
  invisible(m)
}

#' @export
#' @rdname tool
rm_tool <- function(m, tool_pos) {
  invisible(rm_tool_cpp(m, tool_pos))
}

#' @export
#' @rdname tool
rm_tool <- function(m, tool_pos) {
  invisible(rm_tool_cpp(m, tool_pos))
}

