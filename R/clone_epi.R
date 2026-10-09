#' Copy epiworld objects
#'
#' Objects in `epiworldR` (models, tools, viruses) are pointers to C++
#' objects, so the assignment operator (`<-`) only copies the pointer:
#' changing the "copy" changes the original. `clone_epi()` creates an
#' independent copy of the underlying C++ object instead.
#'
#' @param x An object of class `epiworld_model`, `epiworld_tool`, or
#' `epiworld_virus`.
#' @param ... Further arguments passed to methods (currently unused).
#' @details
#' Copies of tools and viruses keep their type and all of their settings
#' (e.g., names, probabilities, distribution functions, and targets). For
#' example, cloning the vaccine of a measles model (retrieved with
#' [get_tool()]) gives another all-or-nothing vaccine, which the generic
#' [tool()] cannot create.
#'
#' Copies of tools and viruses are *new* objects: once added to a model with
#' [add_tool()] or [add_virus()], they get their own id (and, for viruses,
#' their own lineage), so the model records them separately from the
#' original. Use [set_name_tool()] or [set_name_virus()] to tell them apart
#' in the outputs.
#'
#' `clone_epi()` replaces `clone_model()`, which is deprecated.
#' @returns A copy of `x`, with the same class.
#' @examples
#' model <- ModelSIRCONN(
#'   name = "COVID-19", n = 1000, prevalence = 0.01, contact_rate = 5,
#'   transmission_rate = 0.4, recovery_rate = 0.95
#' )
#'
#' vax <- tool(
#'   name = "Vaccine", prevalence = 0.5, as_proportion = TRUE,
#'   susceptibility_reduction = 0.9, transmission_reduction = 0.5,
#'   recovery_enhancer = 0.5, death_reduction = 0.9
#' )
#'
#' # A second vaccine with a lower efficacy
#' vax2 <- clone_epi(vax)
#' set_name_tool(vax2, "Vaccine (lower efficacy)")
#' set_susceptibility_reduction(vax2, 0.7)
#'
#' add_tool(model, vax)
#' add_tool(model, vax2)
#' model
#' @export
clone_epi <- function(x, ...) UseMethod("clone_epi")

#' @export
#' @rdname clone_epi
clone_epi.epiworld_tool <- function(x, ...) {
  structure(
    clone_tool_cpp(x),
    class = class(x)
  )
}

#' @export
#' @rdname clone_epi
clone_epi.epiworld_virus <- function(x, ...) {
  structure(
    clone_virus_cpp(x),
    class = class(x)
  )
}

#' @export
#' @rdname clone_epi
clone_epi.epiworld_model <- function(x, ...) {
  structure(
    clone_model_cpp(x),
    class = class(x)
  )
}

#' @export
clone_epi.default <- function(x, ...) {
  stop(
    "clone_epi() has no method for objects of class '",
    paste(class(x), collapse = "/"), "'.",
    call. = FALSE
  )
}
