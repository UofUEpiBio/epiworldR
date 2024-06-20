stopifnot_agent <- function(x) {
  if (!inherits(x, "epiworld_agent"))
    stop("x must be an object of class epiworld_agent")
}

#' Load agents to a model
#' 
#' These functions provide access to the network of the model. The network is
#' represented by an edgelist. The `agents_smallworld` function generates a
#' small world network with the Watts-Strogatz algorithm. The 
#' `agents_from_edgelist` function loads a network from an edgelist.
#' The `get_network` function returns the edgelist of the network.
#' 
#' @param model Model object of class [epiworld_model].
#' @param source,target Integer vectors describing the source and target of
#' in the edgelist.
#' @param n,size Number of individuals in the population.
#' @param k Number of ties in the small world network.
#' @param d,directed Logical scalar. Whether the graph is directed or not.
#' @param p Probability of rewiring.
#' @export
#' @return
#' - The 'agents_smallworld' function returns a model with the agents 
#' loaded.
#' @examples
#' 
#' # Initializing SIR model with agents_smallworld
#' sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, transmission_rate = 0.9, 
#'                 recovery_rate = 0.1)
#' agents_smallworld(
#'    sir,
#'    n = 1000, 
#'    k = 5,
#'    d = FALSE,
#'    p = .01
#'  )
#' run(sir, ndays = 100, seed = 1912)
#' sir
#' 
#' # We can also retrieve the network
#' net <- get_network(sir)
#' head(net)
#' 
#' # Simulating a bernoulli graph
#' set.seed(333)
#' n <- 1000
#' g <- matrix(runif(n ^ 2) < .01, nrow = n)
#' diag(g) <- FALSE
#' el <- which(g, arr.ind = TRUE) - 1L
#' 
#' 
#' # Generating an empty model
#' sir <- ModelSIR("COVID-19", .01, .8, .3)
#' agents_from_edgelist(
#'   sir,
#'   source = el[,1],
#'   target = el[,2],
#'   size   = n,
#'   directed = TRUE
#' )
#' 
#' # Running the simulation
#' run(sir, 50)
#' 
#' plot(sir)
agents_smallworld <- function(model, n, k, d, p)
  UseMethod("agents_smallworld")

#' @export
agents_smallworld.epiworld_model <- function(model, n, k, d, p) {
  agents_smallworld_cpp(model, n, k, d, p)
  invisible(model)
}

#' @export
#' @return
#' - The `agents_from_edgelist` function returns an empty model of class
#' `epiworld_model`. 
#' @rdname agents_smallworld
agents_from_edgelist <- function(
  model, source, target, size, directed
  ) UseMethod("agents_from_edgelist")

#' @export
agents_from_edgelist.epiworld_model <- function(
  model, source, target, size, directed
  ) {
  
  agents_from_edgelist_cpp(
    model,
    source,
    target,
    size,
    directed
  )
  
  invisible(model)
  
}

#' @export 
#' @rdname agents_smallworld
#' @aliases network
#' @return
#' - The `get_network` function returns a data frame with two columns
#' (`source` and `target`) describing the edgelist of the network.
get_network <- function(model) {
  stopifnot_model(model)
  get_network_cpp(model)
}

#' @export 
#' @return 
#' - `get_agents_states` returns an character vector with the states of the
#' agents by the end of the simulation.
#' @rdname agents_smallworld
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
#' @rdname agents_smallworld
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
#' @rdname agents_smallworld
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
#' @rdname agents_smallworld
#' @return
#' - The functions `has_virus` and `has_tool` return a logical scalar
#' indicating whether the agent has the virus/tool or not.
has_virus <- function(agent, virus) {
  stopifnot_agent(agent)
  stopifnot_virus(virus)
  has_virus_cpp(agent, virus)
}

#' @export
#' @rdname agents_smallworld
has_tool <- function(agent, tool) {
  stopifnot_agent(agent)
  stopifnot_tool(tool)
  has_tool_cpp(agent, tool)
}

#' @export
#' @rdname agents_smallworld
change_state <- function(agent, model, state_new, queue = -99) {
  stopifnot_agent(agent)
  stopifnot_model(model)
  change_state_cpp(agent, model, state_new, queue)
}
