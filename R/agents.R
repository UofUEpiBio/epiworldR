#' Agent functions
#'
#' Functions to retrieve and modify the agents of a model.
#'
#' @param model Model object of class [epiworld_model].
#' @export
#' @concept model-utility-functions
#' @name agent-functions
#' @return
#' - `get_agents_states` returns an character vector with the states of the
#' agents by the end of the simulation.
#' @examples
#' sir <- ModelSIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_rate = 0.9, recovery_rate = 0.1)
#' add_layer(sir, layer_network("contacts", graph_smallworld(k = 5, p = .01)))
#' run(sir, ndays = 100, seed = 1912)
#' table(get_agents_states(sir))
get_agents_states <- function(model) {
  stopifnot_model(model)
  get_agents_states_cpp(model)
}

#' @export
#' @param agent Agent object of class `epiworld_agent`.
#' @param virus Virus object of class `epiworld_virus`.
#' @param state_new Integer scalar. New state of the agent after the action is executed.
#' @param queue Integer scalar. Change in the queuing system after the action is executed.
#' @details
#' The `new_state` and `queue` parameters are optional. If they are not
#' provided, the agent will be updated with the default values of the virus/tool.
#' @rdname agent-functions
#' @return
#' - The function `add_virus_agent` adds a virus to an agent and
#' returns the agent invisibly.
add_virus_agent <- function(
  agent, model, virus, state_new = -99, queue = -99
) {

  stopifnot_model(model)
  stopifnot_virus(virus)
  stopifnot_agent(agent)

  invisible(
    add_virus_agent_cpp(agent, model, virus, state_new, queue)
  )

}

#' @export
#' @param tool Tool object of class `epiworld_tool`.
#' @rdname agent-functions
#' @return
#' - The function `add_tool_agent` adds a tool to an agent and
#' returns the agent invisibly.
add_tool_agent <- function(
  agent, model, tool, state_new = -99, queue = -99
) {

  stopifnot_model(model)
  stopifnot_tool(tool)
  stopifnot_agent(agent)

  invisible(
    add_tool_agent_cpp(agent, model, tool, state_new, queue)
  )

}

#' @export
#' @rdname agent-functions
#' @return
#' - The functions `has_virus` and `has_tool` return a logical scalar
#' indicating whether the agent has the virus/tool or not.
has_virus <- function(agent, virus) {
  stopifnot_agent(agent)
  stopifnot_virus(virus)
  has_virus_cpp(agent, virus)
}

#' @export
#' @rdname agent-functions
has_tool <- function(agent, tool) {
  stopifnot_agent(agent)
  stopifnot_tool(tool)
  has_tool_cpp(agent, tool)
}

#' @export
#' @rdname agent-functions
change_state <- function(agent, model, state_new, queue = -99) {
  stopifnot_agent(agent)
  stopifnot_model(model)
  change_state_cpp(agent, model, state_new, queue)
}
