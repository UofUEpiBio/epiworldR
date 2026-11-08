#' Agent functions
#'
#' Functions to retrieve and modify the agents of a model.
#'
#' @param model Model object of class [epiworld_model].
#' @export
#' @concept model-utility-functions
#' @name agent-functions
#' @return
#' - `get_agents_stages` and `get_agents_statuses` return a character vector
#' with each agent's stage (status), by the end of the simulation.
#' @examples
#' sir <- ModelSIR(name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_prob = 0.9, recovery_rate = 0.1)
#' add_layer(sir, layer_network("contacts", graph_smallworld(k = 5, p = .01)))
#' run(sir, ndays = 100, seed = 1912)
#' table(get_agents_stages(sir))
get_agents_stages <- function(model) {
  stopifnot_model(model)
  get_agents_stages_cpp(model, FALSE)
}

#' @export
#' @rdname agent-functions
get_agents_statuses <- function(model) {
  stopifnot_model(model)
  get_agents_stages_cpp(model, TRUE)
}

#' @export
#' @param agent Agent object of class `epiworld_agent`.
#' @param virus Virus object of class `epiworld_virus`.
#' @details
#' `add_virus_agent()` takes the agent through the virus's infection
#' transition (into the infection's target), as a contact would; the agent must
#' be in a stage the infection leaves. The change applies at the end of the
#' day, with the model's other events. To move an agent along any other
#' declared transition, use [fire()].
#'
#' The virus and the tool are the model's own copies, from [get_virus()] and
#' [get_tool()] (the copies given to [add_virus()] and [add_tool()] are not in
#' the model). These functions are typically called from a global event (see
#' [globalevent_fun()]), with the agent from `get_agents(model)[i]`.
#' @rdname agent-functions
#' @return
#' - The function `add_virus_agent` adds a virus to an agent and
#' returns the agent invisibly.
add_virus_agent <- function(agent, model, virus) {

  stopifnot_model(model)
  stopifnot_virus(virus)
  stopifnot_agent(agent)

  invisible(
    add_virus_agent_cpp(agent, model, virus)
  )

}

#' @export
#' @param tool Tool object of class `epiworld_tool`.
#' @rdname agent-functions
#' @return
#' - The function `add_tool_agent` adds a tool to an agent and
#' returns the agent invisibly.
add_tool_agent <- function(agent, model, tool) {

  stopifnot_model(model)
  stopifnot_tool(tool)
  stopifnot_agent(agent)

  invisible(
    add_tool_agent_cpp(agent, model, tool)
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
