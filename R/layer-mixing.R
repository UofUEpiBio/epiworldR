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
#' `layer_mixing()` creates one to add to any model with [add_layer()]. Every
#' agent of a group takes part in the mixing (the contacts are spread over all
#' of them), and the agents in the model's infectious states transmit (see
#' [set_infectious_states()]; e.g., Infected, not Exposed, in [ModelSEIR()]).
#'
#' With `set_contact_rate()`, the contacts of a layer with a contact matrix are
#' a model parameter times the matrix. The parameter is read at the end of
#' every day, so a global event or a calibration that changes it changes the
#' next day's contacts. With everyone in one group and a 1 x 1 matrix of 1,
#' every agent meets "Contact rate" others a day, at random (what the former
#' connected models, e.g. `ModelSIRCONN()`, did; see the example).
#'
#' @param name String. Name of the layer (unique within a model).
#' @param groups Integer vector with the group of each agent (starting at 0),
#' or `NA` to leave the agent out of the layer. Its length is the population
#' size of the model the layer is added to.
#' @param contact_matrix Numeric square matrix with one row per group.
#' @param contact_prob Numeric vector with the probability that two members of
#' a group meet on a given day: one per group, or one for all groups. Each in
#' `[0, 1]`.
#' @param layer A mixing layer (class `epiworld_layer_mixing`), from
#' `layer_mixing()` or [get_layer()].
#' @param agent Integer. Id of an agent (starting at 0).
#' @param group Integer. Id of a group (starting at 0); in `set_group()`, `NA`
#' takes the agent out of the layer.
#' @param param String. In `set_contact_rate()`, the name of the model
#' parameter the contact matrix is multiplied by (add it with [add_param()]).
#' @param from,to Integers. In `groupmates_of()`, the first and last day of
#' the window (days after the last one simulated are left out).
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
#' A layer with contact probabilities that is traced ([enable_tracing()],
#' before the run) also keeps the run's history of who was in which group and
#' who was absent (in a state the model takes out of the layers, e.g., isolated
#' in [ModelSEIRQuarantine()]). `groupmates_of()` reads it: who shared a group
#' with `agent`, both present, between days `from` and `to`. With
#' `contact_prob = 1` that is exactly who the agent met; below 1, who it could
#' have met ("notify the household"). The history covers whole runs: tracing
#' turned on during a run, or a grouping edited after it, gives none (an
#' error) until the next run.
#'
#' @return
#' - `layer_mixing()` returns a mixing layer, of class
#' `epiworld_layer_mixing`.
#' - `get_contact_prob()` returns the contact probabilities (one per group once
#' the model has run; empty for a layer with a contact matrix).
#' - `set_contact_prob()`, `set_contact_rate()` and `set_group()` return the
#' layer invisibly.
#' - `get_group()` returns the group of the agent (`NA` if none).
#' - `members_of()` returns the ids of the agents in the group, in no
#' particular order.
#' - `groupmates_of()` returns a data frame with one row per groupmate, in
#' ascending id order: `agent`, `last_day` (the last day they shared a group,
#' both present), and `n_days` (how many days they did).
#' @examples
#' # 1,000 agents in households of 4: everyone meets their housemates daily
#' sir <- ModelSIR(
#'   name              = "flu",
#'   n                 = 1000,
#'   prevalence        = .01,
#'   transmission_prob = .2,
#'   recovery_rate     = .1
#' )
#'
#' households <- layer_mixing(
#'   "households",
#'   groups       = (0:999) %/% 4,
#'   contact_prob = 1
#' )
#' add_layer(sir, households)
#'
#' # A daily contact matrix between two groups, children and adults
#' cm <- matrix(c(10, 3, 3, 6), nrow = 2)
#' add_layer(
#'   sir,
#'   layer_mixing(
#'     "community",
#'     groups         = rep(0:1, c(200, 800)),
#'     contact_matrix = cm
#'   )
#' )
#'
#' enable_tracing(households)
#' run(sir, ndays = 50, seed = 1)
#' members_of(households, 0)
#' get_contact_prob(households)
#'
#' # Who shared a household with agent 0 in the last week
#' groupmates_of(households, 0, from = 44, to = 50)
#'
#' # Random mixing: everyone meets "Contact rate" others a day
#' seir <- ModelSEIR(
#'   name              = "flu",
#'   n                 = 1000,
#'   prevalence        = .01,
#'   transmission_prob = .1,
#'   latent_period     = 3,
#'   recovery_rate     = 1 / 7
#' )
#' add_param(seir, "Contact rate", 4)
#'
#' everyone <- layer_mixing("everyone", groups = rep(0, 1000), contact_matrix = 1)
#' set_contact_rate(everyone, "Contact rate")
#' add_layer(seir, everyone)
#'
#' run(seir, ndays = 50, seed = 1)
#' @name mixing-layers
#' @concept contact-layers
NULL

#' @export
#' @rdname mixing-layers
layer_mixing <- function(
  name,
  groups,
  contact_matrix = NULL,
  contact_prob   = NULL
) {

  stopifnot_string(name)

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

  ptr <- layer_mixing_cpp(name, as.integer(groups), as.integer(n_groups))

  if (!is.null(contact_matrix)) {
    stopifany_na(contact_matrix)
    set_contact_matrix_cpp(ptr, as.vector(contact_matrix))
  } else {
    stopifnot_double(contact_prob, lb = 0, ub = 1)
    set_contact_prob_cpp(ptr, as.double(contact_prob))
  }

  new_layer(ptr, "mixing")

}

#' @export
#' @rdname mixing-layers
set_contact_rate <- function(layer, param) {
  stopifnot_layer(layer, "mixing")
  stopifnot_string(param)
  invisible(set_contact_rate_cpp(layer, param))
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
  stopifnot_int(agent, lb = 0)
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
  stopifnot_int(group, lb = 0)
  members_of_cpp(layer, group)
}

#' @export
#' @rdname mixing-layers
groupmates_of <- function(
  layer,
  agent,
  from = 1L,
  to   = .Machine$integer.max
) {

  stopifnot_layer(layer, "mixing")
  stopifnot_int(agent, lb = 0)
  stopifnot_int(from)
  stopifnot_int(to)
  groupmates_of_cpp(layer, agent, from, to)

}
