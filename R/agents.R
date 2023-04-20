#' Load agents to a model
#' @param m Model object.
#' @param source,target Integer vectors describing the source and target of
#' in the edgelist.
#' @param n,size Number of individuals in the population.
#' @param k Number of ties in the small world network.
#' @param d,directed Logical scalar. Whether the graph is directed or not.
#' @param p Probability of rewiring.
#' @export
#' @aliases agents
#' @examples
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
agents_smallworld <- function(m, n, k, d, p)
  UseMethod("agents_smallworld")

#' @export
agents_smallworld.epiworld_model <- function(m, n, k, d, p) {
  agents_smallworld_cpp(m, n, k, d, p)
  invisible(m)
}

#' @export
#' @rdname agents_smallworld
agents_from_edgelist <- function(
  m, source, target, size, directed
  ) UseMethod("agents_from_edgelist")

#' @export
agents_from_edgelist.epiworld_model <- function(
  m, source, target, size, directed
  ) {
  
  agents_from_edgelist_cpp(
    m,
    source,
    target,
    size,
    directed
  )
  
  invisible(m)
  
}

  