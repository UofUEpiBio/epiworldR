#' Mixing layers
#'
#' A mixing layer draws contacts afresh every day, by group, instead of
#' keeping ties. Its contacts are given in one of two ways:
#' - A contact matrix, for a few groups that meet each other (e.g., age groups
#' in public spaces): entry `(g, h)` is the expected number of contacts per day
#' of an agent in group `g` with agents in group `h`.
#' - A contact probability per pair, for many small groups that mix only within
#' (households, classrooms): each member meets each other member of its group
#' with probability `contact_prob` per day, whatever the size of the group.
#' `contact_prob = 1` makes every group a clique (everyone meets everyone,
#' every day), with no ties stored.
#'
#' The mixing models (e.g., [ModelSEIRMixing()]) create their own mixing layer,
#' named `"mixing"`; `layer_mixing()` creates one to add to any model with
#' [add_layer()].
#'
#' @param name String. Name of the layer (unique within a model).
#' @param groups Integer vector with the group of each agent (starting at 0),
#' or `NA` to leave the agent out of the layer. Its length is the population
#' size of the model the layer is added to.
#' @param infectious_states Integer vector. States (starting at 0) whose agents
#' can be met and transmit.
#' @param available_states Integer vector. States (starting at 0) whose agents
#' take part in the mixing (e.g., not hospitalized). With a contact matrix,
#' the contacts per day are spread over the available agents.
#' @param contact_matrix Numeric square matrix with one row per group.
#' @param contact_prob Numeric vector with the probability that two members of
#' a group meet on a given day: one per group, or one for all groups. Each in
#' `[0, 1]`.
#' @param layer A mixing layer (class `epiworld_layer_mixing`), from
#' `layer_mixing()` or [get_layer()].
#' @param agent Integer. Id of an agent (starting at 0).
#' @param group Integer. Id of a group (starting at 0); in `set_group()`, `NA`
#' takes the agent out of the layer.
#'
#' @details
#' Give exactly one of `contact_matrix` and `contact_prob`. The number of
#' groups is the number of rows of `contact_matrix`, or the length of
#' `contact_prob`; with a single `contact_prob`, it is `max(groups) + 1`.
#'
#' As for every layer, an edit made while the model runs (e.g., from a global
#' event) lasts for that run only, and one made between runs is kept. A move
#' made with `set_group()` during a day takes effect at the end of the day.
#' A layer with a contact matrix and one with contact probabilities cannot be
#' swapped during a run.
#'
#' @return
#' - `layer_mixing()` returns a mixing layer, of class
#' `epiworld_layer_mixing`.
#' - `get_contact_prob()` returns the contact probabilities (one per group once
#' the model has run; empty for a layer with a contact matrix).
#' - `set_contact_prob()` and `set_group()` return the layer invisibly.
#' - `get_group()` returns the group of the agent (`NA` if none).
#' - `members_of()` returns the ids of the agents in the group, in no
#' particular order.
#' @examples
#' # 1,000 agents in households of 4: everyone meets their housemates daily
#' sir <- ModelSIR(
#'   name              = "flu",
#'   n                 = 1000,
#'   prevalence        = .01,
#'   transmission_rate = .2,
#'   recovery_rate     = .1
#' )
#'
#' households <- layer_mixing(
#'   "households",
#'   groups            = (0:999) %/% 4,
#'   infectious_states = 1,
#'   available_states  = 0:2,
#'   contact_prob      = 1
#' )
#' add_layer(sir, households)
#'
#' # A daily contact matrix between two groups, children and adults
#' cm <- matrix(c(10, 3, 3, 6), nrow = 2)
#' add_layer(
#'   sir,
#'   layer_mixing(
#'     "community",
#'     groups            = rep(0:1, c(200, 800)),
#'     infectious_states = 1,
#'     available_states  = 0:2,
#'     contact_matrix    = cm
#'   )
#' )
#'
#' run(sir, ndays = 50, seed = 1)
#' members_of(households, 0)
#' get_contact_prob(households)
#' @name mixing-layers
#' @concept contact-layers
NULL

#' @export
#' @rdname mixing-layers
layer_mixing <- function(
  name,
  groups,
  infectious_states,
  available_states,
  contact_matrix = NULL,
  contact_prob   = NULL
) {

  stopifnot_string(name)
  stopifnot_int(infectious_states, lb = 0)
  stopifnot_int(available_states, lb = 0)

  # A fraction would be truncated, and a negative id would be NO_GROUP in C++
  if (
    !(is.numeric(groups) || all(is.na(groups))) ||
      any(groups < 0 | groups != round(groups), na.rm = TRUE)
  )
    stop("-groups- must be non-negative integers (NA for no group).")

  if (is.null(contact_matrix) == is.null(contact_prob))
    stop("Give exactly one of -contact_matrix- and -contact_prob-.")

  n_groups <- if (!is.null(contact_matrix)) {
    nrow(as.matrix(contact_matrix))
  } else if (length(contact_prob) > 1L) {
    length(contact_prob)
  } else {
    max(groups, -1L, na.rm = TRUE) + 1L
  }

  ptr <- layer_mixing_cpp(
    name,
    as.integer(groups),
    as.integer(n_groups),
    as.integer(infectious_states),
    as.integer(available_states)
  )

  if (!is.null(contact_matrix)) {
    stopifany_na(contact_matrix)
    set_contact_matrix_cpp(ptr, as.vector(contact_matrix))
  } else {
    stopifnot_double(contact_prob, lb = 0, ub = 1)
    set_contact_prob_cpp(ptr, as.double(contact_prob))
  }

  # The number of agents, which get_group() checks against
  attr(ptr, "size") <- length(groups)
  new_layer(ptr, "mixing")

}

#' @export
#' @rdname mixing-layers
get_contact_prob <- function(layer) {
  stopifnot_layer(layer, "mixing")
  get_contact_prob_cpp(layer)
}

#' @export
#' @rdname mixing-layers
set_contact_prob <- function(layer, contact_prob) {
  stopifnot_layer(layer, "mixing")
  stopifnot_double(contact_prob, lb = 0, ub = 1)
  invisible(set_contact_prob_cpp(layer, as.double(contact_prob)))
}

#' @export
#' @rdname mixing-layers
get_group <- function(layer, agent) {

  stopifnot_layer(layer, "mixing")

  # A layer of a mixing model gets its groups when the model runs
  n <- attr(layer, "size")
  if (is.null(n)) {
    if (get_n_groups_cpp(layer) == 0L)
      stop("The layer has no groups yet: they are set when the model runs.")
    n <- size(attr(layer, "model"))
  }

  stopifnot_int(agent, lb = 0, ub = n - 1)
  get_group_cpp(layer, agent)

}

#' @export
#' @rdname mixing-layers
set_group <- function(layer, agent, group) {

  stopifnot_layer(layer, "mixing")
  stopifnot_int(agent, lb = 0)
  if (!is.na(group))
    stopifnot_int(group, lb = 0)

  invisible(set_group_cpp(layer, agent, as.integer(group)))

}

#' @export
#' @rdname mixing-layers
members_of <- function(layer, group) {
  stopifnot_layer(layer, "mixing")
  stopifnot_int(group, lb = 0, ub = get_n_groups_cpp(layer) - 1)
  members_of_cpp(layer, group)
}
