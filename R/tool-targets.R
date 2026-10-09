#' Restrict tools to specific viruses
#'
#' By default, a tool acts on every virus in the model. These functions
#' restrict a tool (all four of its effects: susceptibility reduction,
#' transmission reduction, recovery enhancer, and death reduction) to specific
#' viruses, e.g., a vaccine that only protects against one disease.
#'
#' @param tool An object of class [epiworld_tool].
#' @param virus An object of class [epiworld_virus], a list of them, or an
#' integer vector of virus lineage ids (see details). In `get_lineage_virus()`,
#' a single [epiworld_virus].
#' @param viruses Same as `virus`. Use `NULL` (or an empty vector) to make the
#' tool act on every virus again.
#'
#' @details
#' Tools target viruses by *lineage*. A lineage is a virus added to a model
#' with [add_virus()] (or created by the model itself, as in
#' [ModelSIRCONN()]) together with all of its mutations. The lineage id is the
#' virus id assigned when the virus was added to the model, and every variant
#' of the virus keeps it, so a tool targeting a virus also acts on its
#' variants. Only lineages 0 to 62 can be targeted.
#'
#' A virus only has a lineage id after it has been added to a model, so call
#' [add_virus()] before passing the virus to `add_target_tool()` or
#' `set_targets_tool()`. Viruses that the model creates itself can be
#' retrieved with [get_virus()].
#'
#' The model keeps its own copy of each tool, so set the targets *before*
#' calling [add_tool()], or modify the model's copy retrieved with
#' [get_tool()].
#'
#' @returns
#' - `add_target_tool()`, `set_targets_tool()`, and `clear_targets_tool()`
#' return the tool (invisibly).
#' - `get_targets_tool()` returns an integer vector with the targeted lineage
#' ids. An empty vector means the tool acts on every virus.
#' - `get_lineage_virus()` returns the lineage id of the virus (an integer), or
#' `NA` if the virus has not been added to a model yet.
#'
#' @examples
#' # A model with two diseases: "Flu" (created by the model) and "Measles"
#' model <- ModelSIRCONN(
#'   name              = "Flu",
#'   n                 = 2000,
#'   prevalence        = 0.01,
#'   contact_rate      = 4,
#'   transmission_rate = 0.5,
#'   recovery_rate     = 0.2
#' )
#'
#' measles <- virus(
#'   name           = "Measles",
#'   prevalence     = 20,
#'   as_proportion  = FALSE,
#'   prob_infecting = 0.5,
#'   recovery_rate  = 0.2
#' )
#'
#' add_virus(model, measles) # Assigns the lineage id
#' get_lineage_virus(get_virus(model, 0)) # Flu: 0
#' get_lineage_virus(measles) # Measles: 1
#'
#' # A vaccine that only protects against measles
#' mmr <- tool(
#'   name                     = "MMR",
#'   prevalence               = 0.5,
#'   as_proportion            = TRUE,
#'   susceptibility_reduction = 0.97,
#'   transmission_reduction   = 0,
#'   recovery_enhancer        = 0,
#'   death_reduction          = 0
#' )
#'
#' add_target_tool(mmr, measles)
#' get_targets_tool(mmr)
#'
#' add_tool(model, mmr)
#' run(model, ndays = 50, seed = 1912)
#' summary(model)
#'
#' # Targets can also be given as lineage ids
#' set_targets_tool(mmr, c(0, 1))
#' get_targets_tool(mmr)
#'
#' # Back to acting on every virus
#' clear_targets_tool(mmr)
#' get_targets_tool(mmr)
#'
#' @export
#' @concept tool-functions
#' @name tool-targets
#' @aliases add_target_tool
add_target_tool <- function(tool, virus) {

  stopifnot_tool(tool)

  for (id in as_lineage_ids(virus))
    add_target_tool_cpp(tool, id)

  invisible(tool)

}

#' @export
#' @rdname tool-targets
set_targets_tool <- function(tool, viruses) {

  stopifnot_tool(tool)

  ids <- if (is.null(viruses)) integer(0) else as_lineage_ids(viruses)

  invisible(set_targets_tool_cpp(tool, ids))

}

#' @export
#' @rdname tool-targets
get_targets_tool <- function(tool) {

  stopifnot_tool(tool)
  as.integer(get_targets_tool_cpp(tool))

}

#' @export
#' @rdname tool-targets
clear_targets_tool <- function(tool) {

  stopifnot_tool(tool)
  invisible(clear_targets_tool_cpp(tool))

}

#' @export
#' @rdname tool-targets
get_lineage_virus <- function(virus) {

  stopifnot_virus(virus)

  id <- get_lineage_virus_cpp(virus)
  if (id < 0) NA_integer_ else id

}

# Turns viruses (one or a list) or numeric lineage ids into lineage ids
as_lineage_ids <- function(x) {

  if (inherits(x, "epiworld_virus"))
    x <- list(x)

  if (is.list(x)) {

    return(vapply(x, function(v) {

      stopifnot_virus(v)

      id <- get_lineage_virus_cpp(v)
      if (id < 0)
        stop(
          "The virus \"", get_name_virus_cpp(v), "\" has no lineage id. ",
          "Add it to a model with add_virus() before targeting it."
        )

      id

    }, integer(1)))

  }

  if (!is.numeric(x) || anyNA(x) || any(x != round(x)))
    stop(
      "Virus targets must be epiworld_virus objects or integer lineage ids. ",
      "The object passed is of class(es): ", paste(class(x), collapse = ", ")
    )

  if (any(x < 0) || any(x > 62))
    stop(
      "Only virus lineages 0 to 62 can be targeted, but got: ",
      paste(x[x < 0 | x > 62], collapse = ", ")
    )

  as.integer(x)

}
