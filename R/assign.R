#' Placement: who gets what at the start of each run
#'
#' One vocabulary places viruses, tools, initial stages and the groups of
#' mixing layers. A placement gives each eligible agent a **label**, or none,
#' drawn with the model's random numbers at the start of every run.
#'
#' @param p Numeric scalar. A share of the eligible agents, in `[0, 1]`.
#' @param n Numeric scalar. A number of the eligible agents.
#' @param x Numeric vector of shares (or counts), one per label. Its names are
#' the labels (stages); without names, the labels are `0`, `1`, ... (groups).
#' @param layer Character scalar. The name of a mixing layer of the model.
#' @param ids Integer vector. The agents (starting at 0).
#' @param labels Character vector (stages; `NA` for none) or integer vector
#' (groups, starting at 0; `NA` for none), one per agent.
#' @param fun A function of `(agent, model)` returning the agent's label (a
#' stage, a group, or `NA` for none).
#'
#' @details
#' | Form | Labels |
#' |---|---|
#' | `assign_share(p)`, `assign_count(n)` | A random share (or number) of the eligible agents |
#' | `assign_shares(x)`, `assign_counts(x)` | Several labels, each its share of the eligible agents (or number), drawn in order from those not yet labeled |
#' | `assign_group_shares(layer, x)`, `assign_group_counts(layer, x)` | A share (or number) of each group of a mixing layer |
#' | `assign_agents(ids)` | The given agents |
#' | `assign_labels(labels)` | `labels[i + 1]` for agent `i`; fixed (the same every run) |
#' | `assign_by(fun)` | `fun(agent, model)`, once per agent and run |
#'
#' What a label means, and who is eligible, depends on the use:
#'
#' | Use | Label | Eligible |
#' |---|---|---|
#' | [add_virus()] | The stage the agent starts in (a carrier stage); a form without labels uses the infection's target | No virus, and a stage that can be infected |
#' | [add_tool()] | -- | Everyone |
#' | [set_initial()] | A stage without the virus | No virus, in the first stage |
#' | [layer_mixing()] (`groups`) | A group | Everyone |
#'
#' At the start of each run the layers are placed first (so a placement can
#' use their groups), then the viruses, the tools and the initial stages.
#'
#' `assign_by()` calls R once per agent at the start of each run, so
#' [run_multiple()] needs `nthreads = 1`.
#'
#' @returns An object of class `epiworld_assign`.
#' @examples
#' flu <- virus("Flu")
#' add_stage(flu, "Susceptible")
#' add_stage(flu, "Infected", infectious = TRUE)
#' add_stage(flu, "Recovered")
#' add_infection(flu, "Susceptible", "Infected", rate = "Transmission rate")
#' add_rate(flu, "Infected", "Recovered", rate = "Recovery rate")
#'
#' model <- Model(1000)
#' add_param(model, "Transmission rate", 0.3)
#' add_param(model, "Recovery rate", 0.1)
#'
#' # 1% of the agents start infected, a fifth start recovered
#' add_virus(model, flu, assign_share(0.01))
#' set_initial(model, assign_shares(c(Recovered = 0.2)))
#'
#' # Agents 0 to 9 get the tool
#' add_tool(model, tool("Vaccine", 0.5, 0, 0, 0), assign_agents(0:9))
#' @name assign
#' @concept assign
NULL

new_assign <- function(form, ...) {
  structure(list(form = form, ...), class = "epiworld_assign")
}

# Objects that call R from C++ (value_fun(), assign_by()) make a model
# single-threaded: run_multiple() then uses one thread (see make_saver.R)
uses_r <- function(...) {
  any(vapply(list(...), function(x) isTRUE(attr(x, "single_threaded")), TRUE))
}

stopifnot_shares <- function(x, what) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(x < 0))
    stop("-", what, "- must be non-negative numbers.", call. = FALSE)
}

# The labels of a vector of shares or counts: its names (stages), or 0, 1, ...
assign_labels_of <- function(x) {
  if (is.null(names(x)))
    return(seq_along(x) - 1L)
  if (any(names(x) == ""))
    stop("Name every share (or count), or none.", call. = FALSE)
  names(x)
}

#' @export
#' @rdname assign
assign_share <- function(p) {
  stopifnot_double(p)
  if (p < 0 || p > 1)
    stop("-p- must be in [0, 1].", call. = FALSE)
  new_assign("share", values = as.double(p))
}

#' @export
#' @rdname assign
assign_count <- function(n) {
  stopifnot_double(n)
  if (n < 0 || n != round(n))
    stop("-n- must be a whole, non-negative number.", call. = FALSE)
  new_assign("count", values = as.double(n))
}

#' @export
#' @rdname assign
assign_shares <- function(x) {
  stopifnot_shares(x, "x")
  if (any(x > 1))
    stop("-x- must be shares in [0, 1].", call. = FALSE)
  new_assign("shares", values = as.double(x), labels = assign_labels_of(x))
}

#' @export
#' @rdname assign
assign_counts <- function(x) {
  stopifnot_shares(x, "x")
  if (any(x != round(x)))
    stop("-x- must be whole numbers.", call. = FALSE)
  new_assign("counts", values = as.double(x), labels = assign_labels_of(x))
}

#' @export
#' @rdname assign
assign_group_shares <- function(layer, x) {
  stopifnot_string(layer)
  stopifnot_shares(x, "x")
  new_assign("group_shares", layer = layer, values = as.double(x))
}

#' @export
#' @rdname assign
assign_group_counts <- function(layer, x) {
  stopifnot_string(layer)
  stopifnot_shares(x, "x")
  if (any(x != round(x)))
    stop("-x- must be whole numbers.", call. = FALSE)
  new_assign("group_counts", layer = layer, values = as.double(x))
}

#' @export
#' @rdname assign
assign_agents <- function(ids) {
  if (!is.numeric(ids) || anyNA(ids) || any(ids < 0 | ids != round(ids)))
    stop("-ids- must be agent ids (whole numbers, starting at 0).", call. = FALSE)
  new_assign("agents", ids = as.integer(ids))
}

#' @export
#' @rdname assign
assign_labels <- function(labels) {
  if (is.character(labels)) {
    return(new_assign("labels", labels = labels))
  }
  if (!(is.numeric(labels) || all(is.na(labels))) ||
    any(labels < 0 | labels != round(labels), na.rm = TRUE)) {
    stop(
      "-labels- must be stages (character) or groups (whole numbers ",
      "starting at 0), with NA for none.",
      call. = FALSE
    )
  }
  new_assign("labels", labels = as.integer(labels))
}

#' @export
#' @rdname assign
assign_by <- function(fun) {
  if (!is.function(fun))
    stop("-fun- must be a function of (agent, model).", call. = FALSE)
  structure(new_assign("by", fun = fun), single_threaded = TRUE)
}

#' @export
print.epiworld_assign <- function(x, ...) {
  what <- switch(
    x$form,
    share        = sprintf("a share (%s) of the eligible agents", x$values),
    count        = sprintf("%s of the eligible agents", x$values),
    shares       = paste0(
      "shares by label: ", paste(x$labels, x$values, sep = " = ", collapse = ", ")
    ),
    counts       = paste0(
      "counts by label: ", paste(x$labels, x$values, sep = " = ", collapse = ", ")
    ),
    group_shares = sprintf("a share of each group of the layer \"%s\"", x$layer),
    group_counts = sprintf("a number of each group of the layer \"%s\"", x$layer),
    agents       = sprintf("%d given agents", length(x$ids)),
    labels       = sprintf("one label per agent (%d)", length(x$labels)),
    by           = "a label from a function of (agent, model)"
  )
  cat("Placement:", what, "\n")
  invisible(x)
}

stopifnot_assign <- function(a, null_ok = TRUE) {
  if (null_ok && is.null(a))
    return(invisible(TRUE))
  if (!inherits(a, "epiworld_assign"))
    stop(
      "A placement is made by one of the assign_*() functions (see ",
      "help(\"assign\")).",
      call. = FALSE
    )
  invisible(TRUE)
}
