#' Contact layers
#'
#' Who meets whom in a model comes from its contact layers. A network layer
#' holds ties between agents; a mixing layer (created by the mixing models,
#' e.g., [ModelSIRMixing()]) samples contacts between groups with a contact
#' matrix. A model can have several layers, and susceptible agents are exposed
#' to the infectious agents of every one of them. Network models
#' (e.g., [ModelSIR()]) start with no layers: add one with `add_layer()`.
#'
#' @param name String. Name of the layer (unique within a model).
#' @param graph Either a network generator (`graph_smallworld()`,
#' `graph_bernoulli()`, `graph_sbm()`), which draws a new network at the start
#' of every run and replicate, or a fixed network (`graph_edgelist()`).
#' @param model Model object of class [epiworld_model].
#' @param layer Layer object of class `epiworld_layer`.
#' @param k Integer. Number of ties per agent in the ring lattice of the small
#' world network.
#' @param p Double. In `graph_smallworld()`, the probability of rewiring a tie;
#' in `graph_bernoulli()`, the probability of each tie.
#' @param directed Logical scalar. Whether the network is directed. In a
#' directed network, a tie `source -> target` is kept by the source only:
#' `target` is one of `source`'s contacts, but not the reverse. Since agents
#' are exposed to their own contacts, `target` can infect `source`, not the
#' other way around. To model "`i` infects `j`", give the tie as `j -> i`.
#' @param block_sizes Integer vector with the size of each block of the
#' stochastic block model (SBM). They must add up to the population size.
#' @param mixing_matrix Numeric matrix with the probability of a tie between
#' (and within) blocks of the SBM.
#' @param source,target Integer vectors with the ends of each tie (agent ids
#' start at 0).
#' @param size Integer. Number of agents in the network (must match the
#' population size of the model it is added to).
#' @param param String. Name of the model parameter holding the proportion of
#' ties rewired at the end of each step.
#' @param rewiring A rewiring function from `rewire_degseq()`.
#' @param i,j Integer. Ids of the agents at the ends of the tie (starting at 0).
#'
#' @details
#' A layer object shares its layer with the model it is added to, so editing
#' it (e.g., with `add_edge()`) edits the model's layer. An edit made while
#' the model runs lasts for that run only; an edit made between runs is kept.
#'
#' @return
#' - `layer_network()` returns a network layer, of class
#' `epiworld_layer_network`.
#' - `add_layer()` returns the model invisibly.
#' - `get_layer()` returns the layer of the model with that name.
#' - The `graph_*()` functions return a network generator (class
#' `epiworld_generator`), and `graph_edgelist()` a fixed network (class
#' `epiworld_graph`).
#' - `rewire_degseq()` returns a rewiring function (class
#' `epiworld_rewiring`), which `set_rewiring()` gives to a network layer
#' (returned invisibly).
#' - `add_edge()` and `rm_edge()` return (invisibly) whether the network
#' changed; `has_edge()` whether the tie exists.
#' - `get_edgelist()` returns a data frame with the ties of a network layer
#' (columns `from` and `to`).
#' @examples
#' sir <- ModelSIR(
#'   name = "COVID-19", n = 1000, prevalence = 0.01,
#'   transmission_rate = 0.9, recovery_rate = 0.1
#' )
#'
#' # A small-world network, drawn again at every run
#' contacts <- layer_network("contacts", graph_smallworld(k = 5, p = .01))
#' add_layer(sir, contacts)
#'
#' run(sir, ndays = 100, seed = 1912)
#' head(get_edgelist(contacts))
#'
#' # A fixed network, from an edgelist
#' set.seed(333)
#' g <- matrix(runif(1000^2) < .01, nrow = 1000)
#' diag(g) <- FALSE
#' el <- which(g, arr.ind = TRUE) - 1L
#'
#' sir2 <- ModelSIR("COVID-19", 1000, .01, .8, .3)
#' add_layer(
#'   sir2,
#'   layer_network("contacts", graph_edgelist(el[, 1], el[, 2], 1000, TRUE))
#' )
#' run(sir2, 50)
#' @name contact-layers
#' @concept contact-layers
NULL

new_layer <- function(ptr, kind) {
  structure(ptr, class = c(paste0("epiworld_layer_", kind), "epiworld_layer"))
}

stopifnot_layer <- function(layer, kind = NULL) {

  if (inherits(layer, "epiworld_model"))
    stop(
      "This function takes a layer, not a model. Retrieve the layer with ",
      "get_layer(model, name).",
      call. = FALSE
    )

  cls <- if (is.null(kind)) "epiworld_layer" else
    paste0("epiworld_layer_", kind)

  if (!inherits(layer, cls))
    stop(
      "The -layer- object must be of class '", cls, "'. The object passed ",
      "to the function is of class(es): ", paste(class(layer), collapse = ", "),
      call. = FALSE
    )

}

#' @export
#' @rdname contact-layers
layer_network <- function(name, graph) {

  stopifnot_string(name)

  ptr <- if (inherits(graph, "epiworld_generator")) {
    layer_network_generator_cpp(name, graph)
  } else if (inherits(graph, "epiworld_graph")) {
    layer_network_graph_cpp(name, graph)
  } else {
    stop(
      "-graph- must be a network generator (e.g., graph_smallworld()) or a ",
      "fixed network (graph_edgelist())."
    )
  }

  new_layer(ptr, "network")

}

#' @export
#' @rdname contact-layers
add_layer <- function(model, layer) {
  stopifnot_model(model)
  stopifnot_layer(layer)
  add_layer_cpp(model, layer)
  invisible(model)
}

#' @export
#' @rdname contact-layers
get_layer <- function(model, name) {

  stopifnot_model(model)
  stopifnot_string(name)

  # The layer belongs to the model, which the attribute keeps alive
  ptr <- get_layer_cpp(model, name)
  attr(ptr, "model") <- model
  new_layer(ptr, get_layer_kind_cpp(ptr))

}

#' @export
print.epiworld_layer <- function(x, ...) {
  cat(sprintf(
    "Contact layer '%s' (%s)\n",
    get_layer_name_cpp(x),
    sub("epiworld_layer_", "", class(x)[1L])
  ))
  invisible(x)
}

#' @export
#' @rdname contact-layers
graph_smallworld <- function(k, p, directed = FALSE) {
  stopifnot_int(k)
  stopifnot_double(p, lb = 0, ub = 1)
  stopifnot_bool(directed)
  structure(
    graph_smallworld_cpp(k, p, directed),
    class = "epiworld_generator"
  )
}

#' @export
#' @rdname contact-layers
graph_bernoulli <- function(p, directed = FALSE) {
  stopifnot_double(p, lb = 0, ub = 1)
  stopifnot_bool(directed)
  structure(graph_bernoulli_cpp(p, directed), class = "epiworld_generator")
}

#' @export
#' @rdname contact-layers
graph_sbm <- function(block_sizes, mixing_matrix) {

  stopifnot_numvector(block_sizes)
  stopifany_na(mixing_matrix)
  if (any(block_sizes < 0))
    stop("-block_sizes- must be non-negative.")

  structure(
    graph_sbm_cpp(
      as.integer(block_sizes),
      as.vector(mixing_matrix),
      row_major = FALSE
    ),
    class = "epiworld_generator"
  )

}

#' @export
#' @rdname contact-layers
graph_edgelist <- function(source, target, size, directed = FALSE) {
  stopifnot_int(size)
  stopifnot_bool(directed)
  structure(
    graph_edgelist_cpp(
      as.integer(source),
      as.integer(target),
      size,
      directed
    ),
    class = "epiworld_graph"
  )
}

#' @export
#' @rdname contact-layers
rewire_degseq <- function(param) {
  stopifnot_string(param)
  structure(rewire_degseq_cpp(param), class = "epiworld_rewiring")
}

#' @export
#' @rdname contact-layers
set_rewiring <- function(layer, rewiring) {
  stopifnot_layer(layer, "network")
  if (!inherits(rewiring, "epiworld_rewiring"))
    stop("-rewiring- must come from rewire_degseq().")
  set_rewiring_cpp(layer, rewiring)
  invisible(layer)
}

#' @export
#' @rdname contact-layers
add_edge <- function(layer, i, j) {
  stopifnot_layer(layer, "network")
  invisible(add_edge_cpp(layer, i, j))
}

#' @export
#' @rdname contact-layers
rm_edge <- function(layer, i, j) {
  stopifnot_layer(layer, "network")
  invisible(rm_edge_cpp(layer, i, j))
}

#' @export
#' @rdname contact-layers
has_edge <- function(layer, i, j) {
  stopifnot_layer(layer, "network")
  has_edge_cpp(layer, i, j)
}

#' @export
#' @rdname contact-layers
get_edgelist <- function(layer) {
  stopifnot_layer(layer, "network")
  get_edgelist_cpp(layer)
}
