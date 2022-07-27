#' Epi Functions
#' @param x This is a number
#' @param y This is another number
#' @param z This is the last number 
#' @return The sum of \code{x}, \code{y}, and \code{z}.
#' @examples 
#' first_function(1, 1, 3)
#' first_function(13, 2, 400)

first_function <- function(x,y,z) {
  x + y + z
}

#' @rdname first_function
second_function <- function(x) {
  
}

#' @details
#' The only function you're likely to need from roxygen2 is [roxygenize()]. 
#' Otherwise refer to the vignettes to see how to format the documentation.
#' @keywords internal
"_PACKAGE"

