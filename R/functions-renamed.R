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

#' @param n Deprecated. Either set the prevalence during the tool/virus
#' initialization or use `set_prevalence_tool`/`set_prevalence_virus`.
#' @export 
#' @rdname epiworldR-deprecated
add_tool_n <- function(model, tool, n) {

  .Deprecated(new = "add_tool")

  set_prevalence_tool(tool, n, as_proportion = FALSE)
  add_tool(model, tool)

}

#' @export 
#' @rdname epiworldR-deprecated
add_virus_n <- function(model, virus, n) {

  .Deprecated(new = "add_virus")
  set_prevalence_virus(virus, n, as_proportion = FALSE)
  add_virus(model, virus)

}