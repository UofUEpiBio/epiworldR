#' Restrict tools to specific diseases
#'
#' By default, a tool acts on every disease in the model. These functions
#' restrict a tool (all four of its effects: susceptibility reduction,
#' transmission reduction, recovery enhancer, and death reduction) to some
#' diseases, by name, e.g., a vaccine that only protects against one disease.
#'
#' @param tool An object of class [epiworld_tool].
#' @param disease Character scalar. A disease, by name (see [get_diseases()]).
#' @param diseases Character vector of diseases. Use `NULL` (or an empty vector)
#' to make the tool act on every disease again.
#'
#' @details
#' A disease is a natural history: a virus built with [virus()] (or created by
#' a model, as in [ModelSIR()]) is a disease, and its copies ([strain()],
#' [clone_epi()], and its mutations) are strains of that disease. A tool that
#' targets a disease acts on all of its strains.
#'
#' The names are checked when the model runs (a name the model does not have
#' is an error). **With several diseases, a tool must name the diseases it
#' acts on**; with one, a tool without targets acts on it.
#'
#' The model keeps its own copy of each tool, so set the targets *before*
#' calling [add_tool()], or modify the model's copy retrieved with
#' [get_tool()].
#'
#' @returns
#' - `add_target_tool()`, `set_targets_tool()`, and `clear_targets_tool()`
#' return the tool (invisibly).
#' - `get_targets_tool()` returns a character vector with the targeted
#' diseases. An empty vector means the tool acts on every disease.
#'
#' @examples
#' # A model with two diseases: "Flu" (created by the model) and "RSV"
#' model <- ModelSIR(
#'   name              = "Flu",
#'   n                 = 2000,
#'   prevalence        = 0.01,
#'   transmission_prob = 0.5,
#'   recovery_rate     = 0.2
#' )
#' add_param(model, "Contact rate", 4)
#' add_layer(
#'   model,
#'   layer_mixing("everyone", groups = rep(0, 2000), contact_matrix = 1) |>
#'     set_contact_rate("Contact rate")
#' )
#'
#' rsv <- virus("RSV") |>
#'   add_stage("Susceptible") |>
#'   add_stage("Infected", infectious = TRUE) |>
#'   add_stage("Recovered") |>
#'   add_infection("Susceptible", "Infected", rate = 0.3) |>
#'   add_rate("Infected", "Recovered", rate = 0.1)
#' add_virus(model, rsv, assign_count(20))
#' get_diseases(model)
#'
#' # A vaccine that only protects against RSV
#' vax <- tool(
#'   name                     = "RSV vaccine",
#'   susceptibility_reduction = 0.9,
#'   transmission_reduction   = 0,
#'   recovery_enhancer        = 0,
#'   death_reduction          = 0
#' )
#'
#' add_target_tool(vax, "RSV")
#' get_targets_tool(vax)
#'
#' add_tool(model, vax, assign_share(0.5))
#' run(model, ndays = 50, seed = 1912)
#' summary(model)
#'
#' @export
#' @concept tool-functions
#' @name tool-targets
#' @aliases add_target_tool
add_target_tool <- function(tool, disease) {

  stopifnot_tool(tool)
  stopifnot_disease_name(disease)
  invisible(add_target_tool_cpp(tool, disease))

}

#' @export
#' @rdname tool-targets
set_targets_tool <- function(tool, diseases) {

  stopifnot_tool(tool)

  if (is.null(diseases))
    diseases <- character(0)

  if (!is.character(diseases) || anyNA(diseases))
    stop(
      "-diseases- must be a character vector of disease names (see ",
      "get_diseases()), or NULL. Tools target diseases by name, not viruses or ",
      "lineage ids.",
      call. = FALSE
    )

  invisible(set_targets_tool_cpp(tool, diseases))

}

#' @export
#' @rdname tool-targets
get_targets_tool <- function(tool) {

  stopifnot_tool(tool)
  get_targets_tool_cpp(tool)

}

#' @export
#' @rdname tool-targets
clear_targets_tool <- function(tool) {

  stopifnot_tool(tool)
  invisible(clear_targets_tool_cpp(tool))

}

#' Restrict tool contact effects to specific layers
#'
#' By default, a tool's susceptibility and transmission reductions apply to
#' contacts through every layer. Layer targets restrict these two effects to
#' contacts through the named layers, for example a mask used only at work.
#' Recovery and death effects are unaffected by layer targets. An infection
#' given by code, rather than through a layer, sees only tools without layer
#' targets for its contact effects. Disease targets ([add_target_tool()])
#' still apply.
#'
#' @param tool An object of class [epiworld_tool].
#' @param layer Character scalar. A layer name, not a layer object or id.
#' @param layers Character vector of layer names. Duplicates are dropped.
#' Use `NULL` or an empty character vector to apply contact effects to every
#' layer again.
#'
#' @details
#' Layer names are checked when the model runs. A layer the model does not
#' have is an error. The model keeps its own copy of each tool, so set targets
#' before calling [add_tool()], or modify the model's copy with [get_tool()].
#'
#' @returns
#' - `add_layer_target_tool()`, `set_layer_targets_tool()`, and
#' `clear_layer_targets_tool()` return the tool invisibly.
#' - `get_layer_targets_tool()` returns a character vector of layer names.
#' An empty vector means contact effects apply to every layer.
#'
#' @examples
#' mask <- tool("Mask", 1, 0, 0, 0)
#' add_layer_target_tool(mask, "work")
#' get_layer_targets_tool(mask)
#' clear_layer_targets_tool(mask)
#'
#' @export
#' @concept tool-functions
#' @name tool-layer-targets
add_layer_target_tool <- function(tool, layer) {
  stopifnot_tool(tool)
  if (!is.character(layer) || length(layer) != 1L || is.na(layer))
    stop("-layer- must be a layer name (a character scalar).", call. = FALSE)
  invisible(add_layer_target_tool_cpp(tool, layer))
}

#' @export
#' @rdname tool-layer-targets
set_layer_targets_tool <- function(tool, layers) {
  stopifnot_tool(tool)

  if (is.null(layers))
    layers <- character(0)

  if (!is.character(layers) || anyNA(layers))
    stop(
      "-layers- must be a character vector of layer names, or NULL.",
      call. = FALSE
    )

  invisible(set_layer_targets_tool_cpp(tool, layers))
}

#' @export
#' @rdname tool-layer-targets
get_layer_targets_tool <- function(tool) {
  stopifnot_tool(tool)
  get_layer_targets_tool_cpp(tool)
}

#' @export
#' @rdname tool-layer-targets
clear_layer_targets_tool <- function(tool) {
  stopifnot_tool(tool)
  invisible(clear_layer_targets_tool_cpp(tool))
}
