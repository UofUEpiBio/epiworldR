#' Social bubbles intervention for network models
#'
#' Restricts agents' *effective* contacts to small, disjoint "bubbles" during an
#' outbreak. It is designed for network-based models (e.g. [ModelSEIR] and
#' [ModelSIR] built on an explicit contact network). Rather than editing the
#' contact network, `bubbles()` attaches a tool to every agent that blocks
#' transmission between agents in different bubbles and, optionally, reduces
#' transmission within a bubble to represent social distancing. Because contact
#' weights are uniform in the network models, this is equivalent to removing
#' out-of-bubble contacts.
#'
#' Households are identified by `household_id` (one entry per agent, in agent
#' order). Two flavors form the bubbles:
#'
#' - `"household"`: households are shuffled and grouped into bubbles of
#'   `group_size` households each. `group_size = 1` reproduces a strict,
#'   household-only lockdown.
#' - `"peer"`: each agent selects up to `group_size` peers among its existing
#'   external contacts (a neighbor in a different household); each selection
#'   merges the two households, and bubbles are the resulting connected groups.
#'
#' @param model An object of class [epiworld_model] whose agents/network are
#' already built.
#' @param household_id Integer vector of household labels, one per agent, in
#' agent order (its length must equal the number of agents).
#' @param flavor Character scalar: `"household"` or `"peer"` (see details).
#' @param group_size Integer scalar. Households per bubble (`"household"`) or the
#' maximum number of external peers per agent (`"peer"`).
#' @param transmission_factor Numeric scalar in \[0, 1\]. Within-bubble
#' transmission multiplier (`1` = no distancing, `0.5` = halved, `0` = fully
#' blocked). Cross-bubble transmission is always blocked.
#' @param start_day Integer scalar. First day the policy is in effect.
#' @param end_day Integer scalar. Day the policy ends (exclusive); negative means
#' it never ends.
#' @param rewire_every Integer scalar. Re-randomize the bubbles every this-many
#' days (for "changing" policies); `0` keeps fixed bubbles.
#' @param name Character scalar. Name used for the intervention's tool and event.
#'
#' @details
#' The bubble partition is (re)computed by a daily scheduler event and takes
#' effect the following simulation step. The intervention uses a single shared
#' state, so when running replicates with [run_multiple] you must use a single
#' thread (`nthreads = 1`).
#'
#' The travel component and per-agent sport partners described in the modelling
#' protocol are not yet implemented.
#'
#' @returns Invisibly returns the (modified in place) `model`.
#'
#' @export
#' @concept global-events
#' @seealso [global-events], [ModelSEIR]
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
#' # Households of ~3 agents (in agent order)
#' household_id <- rep(seq_len(2000 / 3 + 1), each = 3)[1:2000]
#'
#' # Strategy 1 (household bubbles of 2 households), fixed, from day 10
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

  if (identical(flavor, "household") && group_size < 1)
    stop("For the `household` flavor, `group_size` must be >= 1.")

  bubbles_cpp(
    model,
    household_id,
    flavor,
    as.integer(group_size),
    as.numeric(transmission_factor),
    as.integer(start_day),
    as.integer(end_day),
    as.integer(rewire_every),
    name
  )

  # A single shared state backs the intervention, so parallel replicates must
  # run single-threaded.
  attr(model, "single_threaded") <- TRUE

  invisible(model)

}
