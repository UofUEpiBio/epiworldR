#' Agents in epiworldR
#' @param model An object of class [epiworld_model].
#' @param x An object of class [epiworld_agents].
#' @seealso agents
#' @export
#' @aliases epiworld_agents
#' @return The `get_agents` function returns an object of class [epiworld_agents].
#' @examples
#'  
#' model_sirconn <- ModelSIRCONN(
#' name                = "COVID-19",
#' n                   = 10000,
#' prevalence          = 0.01,
#' contact_rate        = 5,
#' prob_transmission   = 0.4,
#' prob_recovery       = 0.95
#' )
#' 
#' run(model_sirconn, ndays = 100, seed = 1912)
#' 
#' x <- get_agents(model_sirconn) # Storing all agent information into object of 
#'                                # class epiworld_agents
#'                              
#' print(x, compressed = FALSE, max_print = 5) # Displaying detailed information of 
#'                                         # the first 5 agents using 
#'                                         # compressed=F. Using compressed=T
#'                                         # results in less-detailed 
#'                                         # information about each agent. 
#'                                         
#' x[0] # Print information about the first agent. Substitute the agent of 
#'      # interest's position where '0' is. 
get_agents <- function(model) {
  
  res <- get_agents_cpp(model)
  
  structure(
    res,
    class = "epiworld_agents",
    model = model
  )
  
}

#' @param x An object of class [epiworld_agents]
#' @param i Index (id) of the agent (from 0 to `n-1`)
#' @export
#' @rdname get_agents
#' @return The `[` method returns an object of class [epiworld_agent].
`[.epiworld_agents` <- function(x, i) {
  
  structure(
    get_agent_cpp(x, i),
    class = "epiworld_agent",
    model = attr(x, "model")
  )
  
}
`[.epiworld_agents` <- function(x, i) {
  
  structure(
    get_agent_cpp(x, i),
    class = "epiworld_agent",
    model = attr(x, "model")
  )
  
}

#' @export
#' @param compressed Logical scalar. When FALSE, it prints detailed information
#' about the agent.
#' @param ... Ignored
#' @rdname get_agents
print.epiworld_agent <- function(x, compressed = FALSE, ...) {
  
  invisible(print_agent_cpp(x, attr(x, "model"), compressed))
  
}

#' @export
#' @param max_print Integer scalar. Maximum number of agents to print.
#' @rdname get_agents
print.epiworld_agents <- function(x, compressed = TRUE, max_print = 10, ...) {
  
  model <- attr(x, "model")
  cat(sprintf("Agents from the model \"%s\":\n", get_name(model)))
  n <- size(model)
  for (i in 1L:min(max_print, n)) {
    
    print(x[i - 1L], compressed)
    
  }
  
  if (n > max_print)
    cat(sprintf("... %i more agents ...\n", n - max_print))
  
  invisible(x)
  
}

#' @export
#' @rdname get_agents
get_state <- function(x) {
  get_state_agent_cpp(x)
}
#'
