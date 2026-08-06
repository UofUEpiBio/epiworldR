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
#' @param max_households Integer scalar, `"peer"` flavor only. The largest number
#' of households a bubble may contain; a nomination that would exceed it is
#' declined. This is what stops one household's choices from chaining into the
#' next (see Details). The default, `2`, represents two households joined by a
#' close contact. Ignored when `flavor = "household"`, where `group_size` is
#' already the cap.
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
#' ## The algorithms
#'
#' Both rules start from the **household contact graph**: one node per
#' household, with an edge between two households whenever at least one member
#' of the first is connected to a member of the second in the agents' contact
#' network. Every random draw uses the model's generator, so a run is
#' reproducible from its seed.
#'
#' `flavor = "household"` grows bubbles from seed households:
#'
#' 1. Visit the households in random order.
#' 2. Skip a household if it already belongs to a bubble. Otherwise open a new
#'    bubble containing it, and take its unassigned neighbours in the household
#'    contact graph as the *frontier*.
#' 3. While the bubble holds fewer than `group_size` households and the frontier
#'    is not empty, draw a household from the frontier at random, add it, and
#'    extend the frontier with that household's unassigned neighbours. A
#'    household tied to several members of the bubble appears in the frontier
#'    more than once and is more likely to be drawn, so stronger ties are
#'    favoured.
#' 4. Stop when no unassigned neighbour is left, even if the bubble is smaller
#'    than `group_size`.
#'
#' Each bubble is therefore a *connected* group of households, though not
#' necessarily one where every pair is directly tied: with `group_size > 2`, two
#' households in the same bubble may be linked only through a third. Households
#' with no available partner stay on their own, so there are always at least
#' `ceiling(n_households / group_size)` bubbles.
#'
#' `flavor = "peer"` has agents choose from the contacts still available:
#'
#' 1. Visit the agents in random order.
#' 2. Skip an agent whose household is already in a full bubble: it has left the
#'    pool and can neither choose nor be chosen.
#' 3. Otherwise draw one of the agent's contacts outside its own household, at
#'    random and without replacement, repeating until the agent has made
#'    `group_size` successful choices or no contact is left that its bubble can
#'    still take in. A draw is accepted when the two households are in different
#'    bubbles *and* the merged bubble would hold at most `max_households`
#'    households; otherwise that contact is unavailable and the agent draws
#'    again. Accepting merges the two households (household commitment: choosing
#'    a peer brings both households into the bubble), and the agent stops as soon
#'    as its bubble is full.
#'
#' The cap is essential. Without it the merges percolate: with a few members per
#' household each choosing someone, the household graph becomes connected and the
#' entire population ends up in a single bubble, imposing no restriction at all.
#'
#' In practice `max_households` is the dial that sets bubble size, while
#' `group_size` rarely binds: since a choice by *any* member commits the whole
#' household, one household's members tend to fill its bubble however many
#' choices each of them is allowed. With `max_households = 2` a single accepted
#' choice fills the bubble, so `group_size` has no effect at all. Households whose
#' every contact was taken first stay on their own.
#'
#' ## Timing and replicates
#'
#' The bubbles are drawn afresh at the start of every run, using that run's seed,
#' and are already in force on day 1. With `rewire_every > 0` they are redrawn
#' during the run, taking effect the following day. Because the intervention
#' keeps a single shared state, replicates with [run_multiple] must use one
#' thread (`nthreads = 1`); `bubbles()` flags the model accordingly.
#'
#' ## What this cannot represent
#'
#' Both rules produce *exclusive* bubbles, as the modelled policies prescribe. A
#' rule that instead gives each person a personal budget of contacts that need be
#' neither mutual nor exclusive -- "up to ten different people a week", say -- is
#' not a partition of the population and cannot be expressed this way.
#'
#' Fixed out-of-bubble "sport partners" and travel between regions are not yet
#' modelled either.
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
#' # One close contact per person, joining the two households:
#' #   bubbles(model, household_id, "peer", group_size = 1)
#' #
#' # As above, but the bubbles are drawn again every week:
#' #   bubbles(model, household_id, "peer", group_size = 1, rewire_every = 7)
#' #
#' # Individually chosen bubbles of up to four households:
#' #   bubbles(model, household_id, "peer", group_size = 2, max_households = 4)
bubbles <- function(
  model,
  household_id,
  flavor              = c("household", "peer"),
  group_size,
  transmission_factor = 1.0,
  start_day           = 0L,
  end_day             = -1L,
  rewire_every        = 0L,
  name                = "Social bubble",
  max_households      = 2L
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
  max_households      <- as_scalar_int(max_households, "max_households")
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

  if (identical(flavor, "peer") && max_households < 2L)
    stop("For the `peer` flavor, `max_households` must be >= 2.")

  bubbles_cpp(
    model,
    household_id,
    flavor,
    group_size,
    transmission_factor,
    start_day,
    end_day,
    rewire_every,
    name,
    max_households
  )

  # A single shared state backs the intervention, so parallel replicates must
  # run single-threaded.
  attr(model, "single_threaded") <- TRUE

  invisible(model)

}
