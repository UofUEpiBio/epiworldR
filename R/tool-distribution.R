#' Tool distribution functions
#'
#' Distribution functions control how a tool is initially assigned across
#' agents in a model.
#'
#' @param tool An object of class `epiworld_tool`.
#' @param distfun An object of class `epiworld_tool_distfun`.
#' @param prevalence Numeric scalar. Prevalence of the tool. In the case of
#' `distribute_tool_to_entities()`, it is a vector of prevalences, one per
#' entity.
#' @param as_proportion Logical scalar. If `TRUE`, `prevalence` is interpreted
#' as a proportion of the total number of agents in the model.
#' @param agents_ids Integer vector. Indices of the agents to which the tool
#' will be assigned.
#' @returns
#' - `set_distribution_tool()` does not return a value. It assigns a
#' distribution function to the specified tool.
#' - `distribute_tool_randomly()` returns a distribution function of class
#' `epiworld_tool_distfun`. When `agents_ids` is not empty, it distributes the
#' tool randomly within that set. Otherwise it uses all agents in the model.
#' - `distribute_tool_to_set()` returns a distribution function of class
#' `epiworld_tool_distfun`.
#' - `distribute_tool_to_entities()` returns a distribution function of class
#' `epiworld_tool_distfun`.
#' @details
#' `set_distribution_tool()` assigns a distribution function to the specified
#' tool of class [epiworld_tool]. Distribution functions can be created with
#' [distribute_tool_randomly()], [distribute_tool_to_set()], and
#' [distribute_tool_to_entities()].
#'
#' `distribute_tool_randomly()` creates a distribution function that randomly
#' assigns the tool to a proportion of the population.
#'
#' `distribute_tool_to_set()` creates a distribution function that assigns the
#' tool to a fixed set of agents.
#'
#' `distribute_tool_to_entities()` creates a distribution function that assigns
#' the tool to a number of agents based on entity-level prevalence. This is
#' only useful for the mixing models.
#'
#' @examples
#' model <- ModelSIRCONN(
#'   name              = "COVID-19",
#'   n                 = 10000,
#'   prevalence        = 0.01,
#'   contact_rate      = 5,
#'   transmission_rate = 0.4,
#'   recovery_rate     = 0.95
#' )
#'
#' vaccine <- tool(
#'   name = "Vaccine",
#'   prevalence = 0.5,
#'   as_proportion = TRUE,
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5,
#'   death_reduction = .9
#' )
#'
#' set_distribution_tool(
#'   vaccine,
#'   distribute_tool_randomly(0.1, TRUE)
#' )
#'
#' add_tool(model, vaccine)
#'
#' set_distribution_tool(
#'   vaccine,
#'   distribute_tool_to_set(1:10)
#' )
#'
#' @name tool-distribution
#' @concept tool-distribution
NULL

#' @export
#' @rdname tool-distribution
set_distribution_tool <- function(tool, distfun) {

  stopifnot_tool(tool)
  stopifnot_tool_distfun(distfun)
  invisible(set_distribution_tool_cpp(tool = tool, distfun = distfun))

}

#' @export
#' @rdname tool-distribution
distribute_tool_randomly <- function(
  prevalence,
  as_proportion,
  agents_ids = integer(0)
) {

  structure(
    distribute_tool_randomly_cpp(
      as.double(prevalence),
      as.logical(as_proportion),
      as.integer(agents_ids)
    ),
    class = "epiworld_tool_distfun"
  )

}

#' @export
#' @rdname tool-distribution
distribute_tool_to_set <- function(
  agents_ids
) {

  structure(
    distribute_tool_to_set_cpp(
      agents_ids
    ),
    class = "epiworld_tool_distfun"
  )

}

#' @export
#' @rdname tool-distribution
distribute_tool_to_entities <- function(
  prevalence,
  as_proportion
) {

  structure(
    distribute_tool_to_entities_cpp(
      as.double(prevalence),
      as.logical(as_proportion)
    ),
    class = "epiworld_tool_distfun"
  )

}
