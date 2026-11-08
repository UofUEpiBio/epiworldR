#' Viruses and their natural history
#'
#' A virus carries its natural history: the **stages** of the disease and the
#' **transitions** between them. Build a virus with `virus()`, declare its
#' stages and transitions, and add it to a model with [add_virus()], which
#' gives the model its stages.
#'
#' @param name Character scalar. The name of the virus.
#' @param virus An object of class `epiworld_virus`.
#' @param label Character scalar. The stage's label.
#' @param infectious Logical scalar. Whether agents in the stage transmit
#' (implies `carrier`).
#' @param carrier Logical scalar. Whether agents in the stage hold the virus.
#' @param stage Character scalar. A stage, by label.
#' @param weight Numeric scalar. The stage's contact weight (default 1). `0`
#' takes its agents out of every layer: they are not met and do not transmit.
#' On networks, a fraction scales the per-contact probability; mixing layers
#' thin the contacts.
#' @param from Character vector. The stages the transition leaves (by label).
#' `add_after()` takes one.
#' @param to Character scalar. The stage the transition enters.
#' @param rate The transition's value: a number, the name of a model parameter,
#' or an object from [value_inverse()], [value_complement()], [value_fun()] or
#' [rate_fun_logit()]. For `add_infection()`, the probability of infection per
#' contact; for `add_rate()`, the daily probability. A number becomes a model
#' parameter (named after the virus and the transition) when the virus is
#' added to a model.
#' @param effect Character scalar. The tool effect a rate takes: `"none"`,
#' `"recovery"` (the probability becomes `1 - (1 - p) * (1 - recovery
#' enhancer)`) or `"death"` (`p * (1 - death reduction)`).
#' @param duration How long the stage lasts, drawn when the agent enters it: an
#' object from one of the [duration] functions.
#' @param since Character scalar or `NULL`. Counts the duration from the day
#' the agent entered this stage (in the same infection), instead of the day it
#' entered `from`.
#' @param transition_name Character scalar or `NULL`. The transition's name
#' (default: `"From -> To"`), used by the outputs and by [fire()].
#' @param param,replacement Character scalars. The model parameter the
#' transitions name, and the one this virus reads instead.
#' @param ... Not used: a virus takes no prevalence or probabilities (see
#' Details).
#'
#' @details
#' Four kinds of transitions leave a stage:
#'
#' - `add_infection()`: by contact, drawn from the model's layers. Taking it
#'   gives the agent the virus. A stage is susceptible because an infection
#'   leaves it.
#' - `add_rate()`: a daily probability. The rates out of a stage compete: one
#'   draw a day picks one of them, or none.
#' - `add_after()`: a duration, drawn when the agent enters the stage; the
#'   agent leaves on the day it ends, unless a rate fired first.
#' - `add_external()`: fired by code, through [fire()].
#'
#' The first stage is where agents start. Entering a stage that is not a
#' carrier releases the virus. A death is a stage with contact weight 0 and no
#' exits.
#'
#' **Strains.** A virus added to a model is a variant with its own id. A strain
#' of the model's virus is a copy (`clone_epi(get_virus(model, 0))`) with its own
#' name; `use_param()` points its transitions to other parameters, so the
#' strain shares the stages and differs in its values.
#'
#' **Placement.** [add_virus()] takes a placement from the [assign] functions:
#' its labels are the stages the agents start in (a carrier stage); a form
#' without labels uses the infection's target.
#'
#' @returns
#' - `virus()` returns a virus, of class `epiworld_virus`.
#' - `add_stage()`, `set_contact()`, `add_infection()`, `add_rate()`,
#'   `add_after()`, `add_external()` and `use_param()` return the virus
#'   (invisibly), so they can be chained.
#' - `get_stages()` returns a data frame with the stages of a virus (label,
#'   whether they carry the virus or transmit it, and their contact weight), or
#'   the labels of a model's stages.
#' @examples
#' # An SIR, from scratch
#' flu <- virus("Flu")
#' add_stage(flu, "Susceptible")
#' add_stage(flu, "Infected", infectious = TRUE)
#' add_stage(flu, "Recovered")
#' add_infection(flu, "Susceptible", "Infected", rate = "Transmission rate")
#' add_rate(
#'   flu, "Infected", "Recovered", rate = "Recovery rate", effect = "recovery"
#' )
#' get_stages(flu)
#'
#' model <- Model(5000)
#' add_param(model, "Transmission rate", 0.1)
#' add_param(model, "Recovery rate", 1 / 7)
#' add_layer(model, layer_network("contacts", graph_smallworld(10, 0.1)))
#' add_virus(model, flu, assign_share(0.01))
#' run(model, ndays = 50, seed = 1)
#' model
#'
#' # A more transmissible strain of the model's virus
#' delta <- clone_epi(get_virus(model, 0))
#' set_name_virus(delta, "Delta")
#' add_param(model, "Delta transmission rate", 0.2)
#' use_param(delta, "Transmission rate", "Delta transmission rate")
#' add_virus(model, delta, assign_count(10))
#' run(model, ndays = 50, seed = 1)
#'
#' # An SEIR with a fixed latent period of 4 days
#' covid <- virus("COVID-19")
#' add_stage(covid, "Susceptible")
#' add_stage(covid, "Exposed", carrier = TRUE)
#' add_stage(covid, "Infected", infectious = TRUE)
#' add_stage(covid, "Recovered")
#' add_infection(covid, "Susceptible", "Exposed", rate = 0.1)
#' add_after(covid, "Exposed", "Infected", duration_fixed(4))
#' add_rate(covid, "Infected", "Recovered", rate = 1 / 7)
#' @export
#' @concept virus-functions
#' @aliases epiworld_virus
virus <- function(name, ...) {

  if (...length())
    stop(
      "virus() takes only a name: declare its stages and transitions ",
      "(add_stage(), add_infection(), add_rate(), add_after()) and place it ",
      "with add_virus(model, virus, assign). See help(\"virus\").",
      call. = FALSE
    )

  stopifnot_string(name)
  structure(virus_cpp(name), class = "epiworld_virus")

}

#' @export
print.epiworld_virus <- function(x, ...) {
  invisible(print_virus_cpp(x))
}

#' @export
#' @rdname virus
add_stage <- function(virus, label, infectious = FALSE, carrier = infectious) {

  stopifnot_virus(virus)
  stopifnot_string(label)
  stopifnot_bool(infectious)
  stopifnot_bool(carrier)
  if (infectious && !carrier)
    stop("An infectious stage is a carrier.", call. = FALSE)

  add_stage_cpp(virus, label, carrier, infectious)
  invisible(virus)

}

#' @export
#' @rdname virus
set_contact <- function(virus, stage, weight) {

  stopifnot_virus(virus)
  stopifnot_string(stage)
  stopifnot_double(weight, lb = 0)
  invisible(set_contact_cpp(virus, stage, as.double(weight)))

}

stopifnot_value <- function(x, what = "rate") {
  ok <- inherits(x, "epiworld_value") ||
    ((is.numeric(x) || is.character(x)) && length(x) == 1L && !is.na(x))
  if (!ok)
    stop(
      "-", what, "- must be a number, the name of a model parameter, or an ",
      "object from value_inverse(), value_complement(), value_fun() or ",
      "rate_fun_logit().",
      call. = FALSE
    )
  if (is.numeric(x))
    x <- as.double(x)
  x
}

stopifnot_transition_name <- function(x) {
  if (!is.null(x))
    stopifnot_string(x)
  invisible(TRUE)
}

#' @export
#' @rdname virus
add_infection <- function(virus, from, to, rate, transition_name = NULL) {

  stopifnot_virus(virus)
  stopifnot_stringvector(from)
  stopifnot_string(to)
  stopifnot_transition_name(transition_name)
  add_infection_cpp(virus, from, to, stopifnot_value(rate), transition_name)
  if (uses_r(rate))
    attr(virus, "single_threaded") <- TRUE
  invisible(virus)

}

#' @export
#' @rdname virus
add_rate <- function(
  virus,
  from,
  to,
  rate,
  effect          = "none",
  transition_name = NULL
) {

  stopifnot_virus(virus)
  stopifnot_stringvector(from)
  stopifnot_string(to)
  stopifnot_string(effect)
  stopifnot_transition_name(transition_name)
  add_rate_cpp(
    virus, from, to, stopifnot_value(rate), effect, transition_name
  )
  if (uses_r(rate))
    attr(virus, "single_threaded") <- TRUE
  invisible(virus)

}

#' @export
#' @rdname virus
add_after <- function(
  virus,
  from,
  to,
  duration,
  since           = NULL,
  transition_name = NULL
) {

  stopifnot_virus(virus)
  stopifnot_string(from)
  stopifnot_string(to)
  if (!inherits(duration, "epiworld_duration"))
    stop(
      "-duration- is made by one of the duration_*() functions.",
      call. = FALSE
    )
  if (!is.null(since))
    stopifnot_string(since)
  stopifnot_transition_name(transition_name)
  add_after_cpp(virus, from, to, duration, since, transition_name)
  if (uses_r(duration))
    attr(virus, "single_threaded") <- TRUE
  invisible(virus)

}

#' @export
#' @rdname virus
add_external <- function(virus, from, to, transition_name) {

  stopifnot_virus(virus)
  stopifnot_stringvector(from)
  stopifnot_string(to)
  stopifnot_string(transition_name)
  add_external_cpp(virus, from, to, transition_name)
  invisible(virus)

}

#' @export
#' @rdname virus
use_param <- function(virus, param, replacement) {

  stopifnot_virus(virus)
  stopifnot_string(param)
  stopifnot_string(replacement)
  invisible(use_param_cpp(virus, param, replacement))

}

#' @export
#' @param x A virus or a model.
#' @rdname virus
get_stages <- function(x) UseMethod("get_stages")

#' @export
get_stages.epiworld_virus <- function(x) {
  as.data.frame(get_stages_virus_cpp(x))
}

#' @export
get_stages.epiworld_model <- function(x) {
  get_stages_model_cpp(x)
}

#' @export
#' @rdname virus
set_name_virus <- function(virus, name) {
  stopifnot_virus(virus)
  invisible(set_name_virus_cpp(virus, name))
}

#' @export
#' @rdname virus
get_name_virus <- function(virus) {
  stopifnot_virus(virus)
  get_name_virus_cpp(virus)
}

# The virus in a model ---------------------------------------------------------

#' @export
#' @rdname virus
#' @param model An object of class `epiworld_model`.
#' @param assign Who gets the virus at the start of each run: a placement from
#' one of the [assign] functions (e.g., `assign_share(0.01)`), or `NULL` for
#' nobody.
add_virus <- function(model, virus, assign = NULL) {

  stopifnot_model(model)
  stopifnot_virus(virus)
  stopifnot_assign(assign)
  add_virus_cpp(model, virus, assign)
  if (uses_r(virus, assign))
    attr(model, "single_threaded") <- TRUE
  invisible(model)

}

#' @export
#' @rdname virus
#' @param virus_pos Integer. The virus's position in the model (starting at 0).
assign_virus <- function(model, virus_pos, assign) {

  stopifnot_model(model)
  stopifnot_assign(assign, null_ok = FALSE)
  assign_virus_cpp(model, as.integer(virus_pos), assign)
  if (uses_r(assign))
    attr(model, "single_threaded") <- TRUE
  invisible(model)

}

#' @export
#' @rdname virus
rm_virus <- function(model, virus_pos) {

  stopifnot_model(model)
  invisible(rm_virus_cpp(model, as.integer(virus_pos)))

}

#' Durations of a stage
#'
#' How many days a stage lasts, drawn when the agent enters it (see
#' [add_after()]). Always a whole number of days, at least one.
#'
#' @param days,mean,sd Each a number, the name of a model parameter, or an
#' object from [value_inverse()] or [value_complement()].
#' @details
#' - `duration_fixed(days)`: `days`, rounded.
#' - `duration_geometric(mean)`: a geometric number of days with that mean (the
#'   same as leaving each day with probability `1 / mean`).
#' - `duration_lognormal(mean, sd)`, `duration_gamma(mean, sd)`: the mean and
#'   standard deviation of the duration, in days. A standard deviation of zero
#'   is fixed.
#' @returns An object of class `epiworld_duration`.
#' @examples
#' duration_fixed(4)
#' duration_lognormal("Incubation period", 2)
#' @name duration
#' @concept virus-functions
NULL

# A duration whose values call R calls R too
new_duration <- function(kind, mean, sd = NULL) {
  res <- duration_cpp(kind, mean, sd)
  if (uses_r(mean, sd))
    attr(res, "single_threaded") <- TRUE
  res
}

#' @export
#' @rdname duration
duration_fixed <- function(days) {
  new_duration("fixed", stopifnot_value(days, "days"))
}

#' @export
#' @rdname duration
duration_geometric <- function(mean) {
  new_duration("geometric", stopifnot_value(mean, "mean"))
}

#' @export
#' @rdname duration
duration_lognormal <- function(mean, sd) {
  new_duration(
    "lognormal", stopifnot_value(mean, "mean"), stopifnot_value(sd, "sd")
  )
}

#' @export
#' @rdname duration
duration_gamma <- function(mean, sd) {
  new_duration(
    "gamma", stopifnot_value(mean, "mean"), stopifnot_value(sd, "sd")
  )
}

#' @export
print.epiworld_duration <- function(x, ...) {
  cat("A stage duration (see help(\"duration\")).\n")
  invisible(x)
}

#' Values of a transition
#'
#' Besides a number or the name of a model parameter, a transition's value (see
#' [add_rate()]) can be one of these.
#'
#' @param param Character scalar. A model parameter.
#' @param fun A function of `(agent, model)` returning the value for that agent.
#' @param vars Integer vector. The columns of the agents' data (starting at 0;
#' see [set_agents_data()]).
#' @param coefs Numeric vector. One coefficient per column.
#' @param model An object of class `epiworld_model` (with the agents' data).
#' @details
#' - `value_inverse(param)`: `1 / param` (e.g., a rate from a mean duration).
#' - `value_complement(param)`: `1 - param`.
#' - `value_fun(fun)`: `fun(agent, model)`. It calls R once per agent and day
#'   it applies to, which is slow, and [run_multiple()] then needs
#'   `nthreads = 1`.
#' - `rate_fun_logit(vars, coefs, model)`: `1 / (1 + exp(-x'b))`, with `x`
#'   the agent's data in the columns `vars` and `b` the coefficients. Computed
#'   in C++.
#' @returns An object of class `epiworld_value`.
#' @examples
#' value_inverse("Infectious period")
#' value_complement("Vaccine efficacy")
#' @name value
#' @concept virus-functions
NULL

#' @export
#' @rdname value
value_inverse <- function(param) {
  stopifnot_string(param)
  value_inverse_cpp(param)
}

#' @export
#' @rdname value
value_complement <- function(param) {
  stopifnot_string(param)
  value_complement_cpp(param)
}

#' @export
#' @rdname value
value_fun <- function(fun) {
  if (!is.function(fun))
    stop("-fun- must be a function of (agent, model).", call. = FALSE)
  res <- value_fun_cpp(fun)
  attr(res, "single_threaded") <- TRUE
  res
}

#' @export
#' @rdname value
rate_fun_logit <- function(vars, coefs, model) {
  stopifnot_model(model)
  rate_fun_logit_cpp(as.integer(vars), as.double(coefs), model)
}

#' @export
print.epiworld_value <- function(x, ...) {
  cat("A transition value (see help(\"value\")).\n")
  invisible(x)
}
