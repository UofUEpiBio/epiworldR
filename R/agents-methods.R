#' Agents in epiworldR
#' @param model An object of class [epiworld_model].
#' @param x An object of class [epiworld_agents].
#' @seealso agents
#' @export
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
#' @rdname get_agents
print.epiworld_agent <- function(x, compressed = FALSE, ...) {
  
  invisible(print_agent_cpp(x, attr(x, "model"), compressed))
  
}

#' @export
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

