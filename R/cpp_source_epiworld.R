#' Extend epiworldR in C++
#'
#' Compile C++ that works on epiworldR's models, in this session. The code
#' reaches a model built in R with `epiworldR::as_model()` and can then use
#' the whole epiworld C++ API on it: add global events, tools, viruses with
#' their natural histories, layers, parameters... See
#' `vignette("extending", package = "epiworldR")`.
#'
#' @param file Path to a C++ file, or `NULL` with `code`.
#' @param code C++ code (a string), or `NULL` with `file`.
#' @param env The environment where the functions marked
#' `[[cpp11::register]]` are defined.
#' @param ... Passed to [cpp11::cpp_source()] (e.g., `quiet`, `clean`).
#' @details
#' The code is compiled with [cpp11::cpp_source()] as C++17, linked to
#' epiworldR's headers, after `#include <epiworldR/extension.h>` (added here,
#' so the code does not need it). That header sets epiworldR's configuration
#' and gives `epiworldR::as_model(SEXP)`, which returns the `epiworld::Model &`
#' behind an R model. What the code adds to a model belongs to it, as if added
#' from R.
#'
#' The epiworld classes must be the same in the code and in the installed
#' epiworldR, so `as_model()` stops with an error if the code was built
#' against other headers or another configuration (rebuild it after updating
#' epiworldR).
#'
#' In a package, add `LinkingTo: cpp11, epiworldR` to `DESCRIPTION`,
#' `CXX_STD = CXX17` to `src/Makevars`, and `#include <epiworldR/extension.h>`
#' to the C++ files.
#' @returns The value of [cpp11::cpp_source()], invisibly: the functions
#' marked `[[cpp11::register]]` are defined in `env`.
#' @export
#' @examples
#' \dontrun{
#' cpp_source_epiworld(code = "
#'   [[cpp11::register]]
#'   int n_agents(SEXP model) {
#'     return static_cast< int >(epiworldR::as_model(model).size());
#'   }
#' ")
#' n_agents(Model(100))
#' }
cpp_source_epiworld <- function(
  file = NULL,
  code = NULL,
  env = parent.frame(),
  ...
) {

  if (!requireNamespace("cpp11", quietly = TRUE)) {
    stop(
      "cpp_source_epiworld() needs the cpp11 package: ",
      "install.packages(\"cpp11\")."
    )
  }

  if (is.null(file) == is.null(code)) {
    stop("Give either `file` or `code`.")
  }

  if (!is.null(file)) {
    code <- paste(readLines(file, warn = FALSE), collapse = "\n")
  }

  if (!is.character(code) || length(code) != 1L) {
    stop("`code` must be a single string.")
  }

  # cpp11 adds epiworldR's headers for the linking_to attribute, which it
  # reads only as code (on a declaration, as [[cpp11::register]] is)
  code <- paste0(
    "[[cpp11::linking_to(\"epiworldR\")]] void epiworldR_extension_();\n",
    "#include <epiworldR/extension.h>\n",
    code
  )

  invisible(cpp11::cpp_source(code = code, env = env, cxx_std = "CXX17", ...))

}
