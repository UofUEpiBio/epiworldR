#' Transmission mode (push or pull)
#'
#' Controls how models compute transmission each day.
#'
#' In a *pull* step, each susceptible agent draws its contacts from every
#' layer: its neighbors in a network, the agents it meets in a mixing layer. In
#' a *push* step, each infected agent adds its infection hazards to the
#' susceptible agents it reaches: its neighbors, and the agents a mixing layer
#' draws for it, each as likely as their own pull would draw it. Both give the
#' same distribution of infections; they differ only in speed and in the
#' random numbers drawn, so the same seed gives different (but equally valid)
#' runs under each mode. A step pushes or pulls every layer.
#'
#' Pushing is cheaper while few agents are infected; pulling can be cheaper
#' near the peak of a large outbreak. `"auto"`, the default, chooses each day:
#' it pushes when \eqn{D_c + 4 N_c \le \kappa (D_s + 4 N_s)}, where \eqn{D} is
#' the sum of the agents' degrees and \eqn{N} the number of agents, over the
#' infected agents that can transmit (\eqn{c}) and the susceptible agents
#' (\eqn{s}), and pulls otherwise. Each mixing layer adds 4 to the price of
#' every agent per group it draws from (all of them with a contact matrix,
#' its own with contact probabilities), on both sides.
#'
#' Only susceptible states that use epiworld's default sampler can be pushed.
#' Models with a directed network always pull. Models with a mixing layer
#' (see [layer_mixing()]) may push too. No mode reproduces runs of earlier
#' versions: since epiworld 1.0.0 every exposure is independent, which
#' changes results in both modes.
#'
#' @param model An `epiworld_model` object.
#' @param mode Character scalar. One of `"auto"`, `"push"`, or `"pull"`.
#' `"pull"` draws every contact from the susceptible side, as epiworld 0.x
#' did.
#' @param kappa Numeric scalar. Threshold used by `"auto"`; a finite,
#' non-negative number. Smaller values pull more often. When `NULL`, the
#' default of the C++ library is used (0.5).
#' @returns
#' - `set_transmission_mode()` returns the model invisibly.
#' - `get_transmission_mode()` returns `"auto"`, `"push"`, or `"pull"`, as set
#'   with `set_transmission_mode()`.
#' - `get_last_transmission_mode()` returns `"push"` or `"pull"`: the mode
#'   used in the most recent step of the last run.
#' - `get_transmission_kappa()` returns the `kappa` threshold.
#' @export
#' @name transmission-mode
#' @examples
#' model <- ModelSIR(
#'   name = "A Virus", 10000, prevalence = .01, transmission_prob = .5,
#'   recovery_rate = 1 / 7
#' )
#' add_layer(model, layer_network("network", graph_smallworld(k = 5, p = .01)))
#' verbose_off(model)
#'
#' get_transmission_mode(model)
#' get_transmission_kappa(model)
#'
#' # Always push
#' set_transmission_mode(model, "push")
#' run(model, ndays = 50, seed = 1912)
#' get_last_transmission_mode(model)
set_transmission_mode <- function(
  model,
  mode  = c("auto", "push", "pull"),
  kappa = NULL
) {

  stopifnot_model(model)
  mode <- match.arg(mode)

  if (is.null(kappa))
    kappa <- default_transmission_kappa_cpp()

  if (length(kappa) != 1L || !is.numeric(kappa) || !is.finite(kappa) ||
    kappa < 0)
    stop("`kappa` must be a single, finite, non-negative number.")

  invisible(set_transmission_mode_cpp(model, mode, as.double(kappa)))

}

#' @export
#' @rdname transmission-mode
get_transmission_mode <- function(model) {
  stopifnot_model(model)
  get_transmission_mode_cpp(model)
}

#' @export
#' @rdname transmission-mode
get_last_transmission_mode <- function(model) {
  stopifnot_model(model)
  get_last_transmission_mode_cpp(model)
}

#' @export
#' @rdname transmission-mode
get_transmission_kappa <- function(model) {
  stopifnot_model(model)
  get_transmission_kappa_cpp(model)
}
