#' Diseases and strains
#'
#' A model can hold several diseases, each with its own natural history (its
#' stages and transitions), on the same population, layers and statuses. An
#' agent carries at most one infection at a time: an agent recovered from one
#' disease can catch another, but an agent with one disease catches nothing
#' else until it leaves the stages that carry the virus.
#'
#' @param model An object of class [epiworld_model].
#' @param virus An object of class [epiworld_virus].
#' @param name Character scalar. The strain's name.
#'
#' @details
#' A virus built with [virus()] (or created by a model, as in [ModelSIR()]) is
#' a disease. A copy of it, from `strain()` (or [clone_epi()], or a mutation),
#' is a strain of the same disease: it shares the stages and transitions, and
#' can read its own parameters with [use_param()]. [add_virus()] adds a strain
#' to its disease, and any other virus as a new disease, which must have a name
#' of its own.
#'
#' With several diseases:
#' - every output (e.g., [get_hist_total()], [get_transmissions()]) returns the
#' rows of every disease, one after the other, with a `disease` column; their
#' `disease` argument selects one;
#' - what refers to a stage needs the disease: [get_stage()],
#' [get_agents_stages()], [get_stages()], [set_initial()] and [fire()] take
#' `disease`, and are an error without it;
#' - tools and layers name the diseases they act on ([add_target_tool()],
#' [add_target_layer()]); a tool without targets is an error;
#' - infections are pulled ([set_transmission_mode()]: `"auto"` pulls, and
#' `"push"` is an error);
#' - when a disease removes agents (a stage with contact weight 0 and no exits,
#' e.g., Dead), the others get a last stage, "Removed (other disease)".
#'
#' @returns
#' - `get_diseases()` returns the names of the model's diseases, in the order
#' they were added (the first is the model's own, in the included models).
#' - `strain()` returns a new virus (of class [epiworld_virus]), a strain of
#' `virus`'s disease named `name`.
#' @examples
#' model <- ModelSIR(
#'   name              = "Flu",
#'   n                 = 2000,
#'   prevalence        = 0.01,
#'   transmission_prob = 0.3,
#'   recovery_rate     = 0.2
#' )
#' add_layer(model, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
#'
#' # A second strain of the flu, more transmissible
#' add_param(model, "Flu B transmission", 0.5)
#' flu_b <- strain(get_virus(model, 0), "Flu B")
#' use_param(flu_b, "Transmission probability", "Flu B transmission")
#' add_virus(model, flu_b, assign_count(10))
#'
#' # A second disease
#' rsv <- virus("RSV") |>
#'   add_stage("Susceptible") |>
#'   add_stage("Infected", infectious = TRUE) |>
#'   add_stage("Recovered") |>
#'   add_infection("Susceptible", "Infected", rate = 0.2) |>
#'   add_rate("Infected", "Recovered", rate = 0.1)
#' add_virus(model, rsv, assign_count(10))
#'
#' get_diseases(model)
#' run(model, ndays = 100, seed = 1)
#' head(get_hist_total(model, disease = "RSV"))
#' table(get_agents_stages(model, disease = "RSV"))
#' @name diseases
#' @concept virus-functions
NULL

#' @export
#' @rdname diseases
get_diseases <- function(model) {
  stopifnot_model(model)
  get_diseases_cpp(model)
}

#' @export
#' @rdname diseases
strain <- function(virus, name) {

  stopifnot_virus(virus)
  stopifnot_string(name)
  structure(
    strain_cpp(virus, name),
    class           = class(virus),
    single_threaded = attr(virus, "single_threaded")
  )

}

# A disease name, as the functions that take one expect
stopifnot_disease_name <- function(disease) {

  if (!is.character(disease) || length(disease) != 1L || is.na(disease))
    stop(
      "-disease- must be a disease's name (a character scalar; see ",
      "get_diseases()).",
      call. = FALSE
    )

  invisible(disease)

}

# The index of a disease (by name) for the C++ functions; NULL is every
# disease (-1)
disease_index <- function(model, disease) {

  if (is.null(disease))
    return(-1L)

  stopifnot_disease_name(disease)

  diseases <- get_diseases(model)
  i <- match(disease, diseases)
  if (is.na(i))
    stop(
      "The model has no disease \"", disease, "\". Its diseases are: ",
      if (length(diseases)) paste0("\"", diseases, "\"", collapse = ", ")
      else "(none)",
      ".",
      call. = FALSE
    )

  i - 1L

}
