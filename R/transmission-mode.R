#' Transmission mode (push or pull)
#'
#' Controls how network models compute transmission each day.
#'
#' In a *pull* step, each susceptible agent looks at its infected neighbors. In
#' a *push* step, each infected agent adds its infection odds to its
#' susceptible neighbors. Both give the same distribution of infections; they
#' differ only in speed and in the random numbers drawn, so the same seed gives
#' different (but equally valid) runs under each mode.
#'
#' Pushing is cheaper while few agents are infected; pulling can be cheaper
#' near the peak of a large outbreak. `"auto"`, the default, chooses each day:
#' it pushes when \eqn{D_c + 4 N_c \le \kappa (D_s + 4 N_s)}, where \eqn{D} is
#' the sum of the agents' degrees and \eqn{N} the number of agents, over the
#' infected agents that can transmit (\eqn{c}) and the susceptible agents
#' (\eqn{s}), and pulls otherwise.
#'
#' Only susceptible states that use epiworld's default network sampler can be
#' pushed. Other models, e.g., mixing models ([ModelSIRMixing()]) and
#' connected models ([ModelSIRCONN()]), always pull, as do directed networks.
#'
#' @param model An `epiworld_model` object.
#' @param mode Character scalar. One of `"auto"`, `"push"`, or `"pull"`.
#' `"pull"` reproduces the random streams of epiworld 0.15 and earlier.
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
#'   name = "A Virus", prevalence = .01, transmission_rate = .5,
#'   recovery_rate = 1/7
#' )
#' agents_smallworld(model, n = 10000, k = 5, d = FALSE, p = .01)
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
