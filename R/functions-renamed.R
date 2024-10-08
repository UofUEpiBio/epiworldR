#' Deprecated and removed functions in epiworldR
#' @description
#' Starting version 0.0-4, epiworld changed how it refered to "actions."
#' Following more traditional ABMs, actions are now called "events."
#'
#' @param ... Arguments to be passed to the new function.
#' @param model Model object of class `epiworld_model`.
#' @param tool Tool object of class `epiworld_tool`.
#' @param virus Virus object of class `epiworld_virus`.
#' @name epiworldR-deprecated
NULL

#' @param n Deprecated.
#' @export
#' @rdname epiworldR-deprecated
add_tool_n <- function(model, tool, n) {

  .Deprecated(new = "add_tool")

  set_distribution_tool(
    tool,
    distfun = distribute_tool_randomly(
      prevalence = n,
      as_proportion = TRUE
    )
  )

  add_tool(model, tool)

}

#' @export
#' @rdname epiworldR-deprecated
add_virus_n <- function(model, virus, n) {

  .Deprecated(new = "add_virus")

  set_distribution_virus(
    virus = virus,
    distfun = distribute_virus_randomly(
      prevalence = n,
      as_proportion = TRUE
    )
  )

  add_virus(model, virus)

}
