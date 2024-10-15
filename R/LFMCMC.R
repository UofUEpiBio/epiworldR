#' LFMCMC
#'
#'
#' @export
LFMCMC <- function() {
  structure(
    LFMCMC_cpp(),
    class = c("epiworld_lfmcmc")
  )
}
