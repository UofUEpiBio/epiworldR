#' Tools in epiworld
#' @param m Model
#' @param name
#' @param susceptibility_reduction
#' @param transmission_reduction
#' @param recovery_enhancer
#' @param death_reduction
#' @examples 
#' epitool <- tool(
#'   "Vaccine",
#'   .9,
#'   .3,
#'   .5,
#'   .9
#' )
#' 
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
#' @rdname tool
add_tool <- function(m, t, prevalence) {
  add_tool_cpp(m, t, prevalence)
  invisible(m)
}

#' @export
#' @rdname tool
add_tool_n <- function(m, t, prevalence) {
  add_tool_n_cpp(m, t, prevalence)
  invisible(m)
}