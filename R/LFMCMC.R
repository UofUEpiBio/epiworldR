#' LFMCMC
#'
#'
#' @export
LFMCMC <- function() {
  structure(
    LFMCMC_cpp(),
    # class = c("epiworld_surv", "epiworld_model")
  )
}
