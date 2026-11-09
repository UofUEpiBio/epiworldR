#' Deprecated and removed functions in epiworldR
#' @description
#' Starting version 0.0-4, epiworld changed how it refered to "actions."
#' Following more traditional ABMs, actions are now called "events."
#'
#' Starting version 0.11.0.0, the measles models have been removed from
#' `epiworldR` and now are part of the `measles` R package:
#' <https://github.com/UofUEpiBio/measles>. The models previously available
#' in `epiworldR` were: `ModelMeaslesSchool` and `ModelMeaslesMixing`.
#'
#' Starting version 0.19.0-0, `clone_model()` is deprecated in favor of
#' [clone_epi()], which also copies tools and viruses.
#'
#' @param ... Arguments to be passed to the new function.
#' @param model Model object of class `epiworld_model`.
#' @name epiworldR-deprecated
#' @concept deprecated-functions
NULL

#' @export
#' @rdname epiworldR-deprecated
clone_model <- function(model) {

  .Deprecated(new = "clone_epi")

  clone_epi(model)

}
