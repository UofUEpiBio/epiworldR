#' epiworldR
#' @useDynLib epiworldR, .registration = TRUE
#' @importFrom graphics boxplot plot
#' @keywords internal
"_PACKAGE"

#' Version of the epiworld C++ code
#'
#' Returns the version of the C++ library epiworld. The code
#' is hosted on GitHub at <https://github.com/UofUEpiBio/epiworld>.
#'
#' @return
#' A character string representing the version of the C++ library.
#'
#' @examples
#' epiworld_cpp_version()
#'
#' @export
#' @keywords internal
epiworld_cpp_version <- function() {
  epiworld_cpp_version_cpp()
}
