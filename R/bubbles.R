#' Social bubbles intervention for network models
#'
#' A "social bubble" policy lets people keep seeing a small, fixed set of others
#' while cutting off the rest of their contacts. `bubbles()` implements such a
#' policy for models built on an explicit contact network (e.g. [ModelSEIR] or
#' [ModelSIR] after [agents_smallworld()] or [agents_from_edgelist()]), following
#' the household-based rules used during COVID-19.
#'
#' @param model An object of class [epiworld_model] whose agents and contact
#' network have already been created.
#' @param household_id Integer vector of household labels, one per agent, in
#' agent order; its length must equal the number of agents. Labels are arbitrary
#' and need not be consecutive. Agents sharing a label form a household and are
#' always placed in the same bubble.
#' @param flavor Character scalar. `"household"` if whole households choose a
#' bubble together, `"peer"` if individuals choose contacts (see Details).
#' @param group_size Integer scalar. For `"household"`, the maximum number of
#' households per bubble (`1` gives a strict household-only lockdown). For
#' `"peer"`, the maximum number of contacts each agent may keep outside its
#' household.
#' @param transmission_factor Numeric scalar in \[0, 1\]. Multiplier applied to
#' transmission *within* a bubble: `1` leaves it untouched, `0.5` halves it
#' (physical distancing), `0` blocks it entirely. Transmission *between* bubbles
#' is always blocked.
#' @param start_day Integer scalar. First day on which the policy applies.
#' @param end_day Integer scalar. Day on which the policy is lifted (exclusive);
#' a negative value means it never lifts. Must be greater than `start_day`.
#' @param rewire_every Integer scalar. Re-randomize the bubbles every this-many
#' days, for policies whose contacts change over time; `0` keeps them fixed.
#' @param name Character scalar. Name given to the intervention's tool and event.
#'
#' @details
#' The contact network is **not** modified. Instead, every agent receives a tool
#' that, on each exposure, compares the bubble of the susceptible agent with the
#' bubble of the infectious one: if they differ, transmission is blocked; if they
#' match, transmission is multiplied by `transmission_factor`. Because contact
#' weights are uniform in these models, this is equivalent to deleting
#' out-of-bubble contacts, while leaving the network available for other purposes
#' (such as contact tracing or network summaries).
#'
#' ## How bubbles are formed
#'
#' Bubbles are always disjoint (each agent is in exactly one), never split a
#' household, and are only ever built between households that are **actually
#' connected** in the contact network.
#'
#' That last property matters: the intervention can only suppress transmission
#' along existing edges, never create new ones, so placing two households that
#' share no contact in the same bubble would change nothing. Pairing households
#' at random would make `group_size` inert -- indistinguishable from a strict
#' lockdown however large the bubbles are. It also matches the real policy: a
#' household picks a bubble partner it already socialises with.
#'
#' The two flavors are:
#'
#' - `"household"`: a bubble grows from one household by repeatedly taking in a
#'   random household connected to it, up to `group_size` households. This is the
#'   household-level rule (e.g. Belgium, 10 May 2020).
#' - `"peer"`: each agent picks up to `group_size` of its existing contacts
#'   outside its household; every pick merges the two households, so a teenager
#'   choosing a partner brings both households into one bubble. This is the
#'   individual-level rule (e.g. Belgium, 19 October 2020, with `group_size = 1`).
#'
#' A household with no available connected partner ends up in a smaller bubble,
#' possibly on its own.
#'
#' ## Timing and replicates
#'
#' The bubbles are drawn afresh at the start of every run, using that run's seed,
#' and are already in force on day 1. With `rewire_every > 0` they are redrawn
#' during the run, taking effect the following day. Because the intervention
#' keeps a single shared state, replicates with [run_multiple] must use one
#' thread (`nthreads = 1`); `bubbles()` flags the model accordingly.
#'
#' Fixed out-of-bubble "sport partners" and travel between regions are not yet
#' modelled.
#'
#' @returns Invisibly returns `model`, which is modified in place.
#'
#' @export
#' @concept global-events
#' @seealso [global-events], [ModelSEIR], [agents_smallworld()]
#' @examples
#' # A network SEIR model on a small-world contact network
#' model <- ModelSEIR(
#'   name              = "Flu",
#'   prevalence        = 0.01,
#'   transmission_rate = 0.1,
#'   incubation_days   = 4.5,
#'   recovery_rate     = 1 / 7
#' )
#'
#' agents_smallworld(model, n = 2000, k = 8, d = FALSE, p = 0.05)
#'
#' # Households of three agents (in agent order)
#' household_id <- rep(seq_len(ceiling(2000 / 3)), each = 3)[1:2000]
#'
#' # Two households per bubble from day 10, halving within-bubble transmission
#' bubbles(
#'   model,
#'   household_id        = household_id,
#'   flavor              = "household",
#'   group_size          = 2,
#'   transmission_factor = 0.5,
#'   start_day           = 10
#' )
#'
#' run(model, ndays = 100, seed = 55)
#' model
#'
#' # Other policies:
#' #
#' # Strict household-only lockdown (a common reference scenario):
#' #   bubbles(model, household_id, "household", group_size = 1)
#' #
#' # One close contact per person, fixed:
#' #   bubbles(model, household_id, "peer", group_size = 1)
#' #
#' # Up to ten contacts per person, renewed weekly:
#' #   bubbles(model, household_id, "peer", group_size = 10, rewire_every = 7)
bubbles <- function(
  model,
  household_id,
  flavor              = c("household", "peer"),
  group_size,
  transmission_factor = 1.0,
  start_day           = 0L,
  end_day             = -1L,
  rewire_every        = 0L,
  name                = "Social bubble"
) {

  stopifnot_model(model)
  flavor <- match.arg(flavor)

  # Coerce/validate the scalar arguments once, up front. Each must be a single,
  # non-missing value; integer-valued arguments must be whole numbers (to avoid
  # silent truncation by as.integer() at the .Call boundary).
  as_scalar_num <- function(x, nm) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x))
      stop(sprintf("`%s` must be a single, non-missing number.", nm))
    as.numeric(x)
  }
  as_scalar_int <- function(x, nm) {
    x <- as_scalar_num(x, nm)
    if (x != round(x))
      stop(sprintf("`%s` must be a whole number.", nm))
    as.integer(x)
  }

  group_size          <- as_scalar_int(group_size, "group_size")
  start_day           <- as_scalar_int(start_day, "start_day")
  end_day             <- as_scalar_int(end_day, "end_day")
  rewire_every        <- as_scalar_int(rewire_every, "rewire_every")
  transmission_factor <- as_scalar_num(transmission_factor, "transmission_factor")

  if (length(name) != 1L || is.na(name) || !is.character(name))
    stop("`name` must be a single, non-missing string.")

  household_id <- as.integer(household_id)
  if (anyNA(household_id))
    stop("`household_id` must contain only non-missing integer household labels.")
  if (any(household_id < 0L))
    stop("`household_id` must be non-negative.")

  n <- size(model)
  if (length(household_id) != n)
    stop(sprintf(
      paste0(
        "`household_id` must have one entry per agent (length %d), ",
        "but its length is %d."
      ),
      n, length(household_id)
    ))

  if (transmission_factor < 0 || transmission_factor > 1)
    stop("`transmission_factor` must be between 0 and 1.")

  if (identical(flavor, "household") && group_size < 1L)
    stop("For the `household` flavor, `group_size` must be >= 1.")

  bubbles_cpp(
    model,
    household_id,
    flavor,
    group_size,
    transmission_factor,
    start_day,
    end_day,
    rewire_every,
    name
  )

  # A single shared state backs the intervention, so parallel replicates must
  # run single-threaded.
  attr(model, "single_threaded") <- TRUE

  invisible(model)

}
