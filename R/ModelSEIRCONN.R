#' Susceptible Exposed Infected Removed model (SEIR connected)
#' @param name Name of the virus
#' @param prevalence a number
#' @param reproductive_number a number
#' @param prob_transmission a number
#' @param incubation_days a number
#' @param prob_recovery a number
#' @export
#' @family Models
ModelSEIRCONN <- function(
    name, n, prevalence, reproductive_number, prob_transmission, incubation_days, prob_recovery
) {
  
  structure(
    ModelSEIRCONN_cpp(name, n, prevalence, reproductive_number, prob_transmission, incubation_days, prob_recovery),
    class = "epiworld_seirconn"
  )
  
}

#' @rdname ModelSEIRCONN
#' @export
init.epiworld_seirconn <- function(m, days, seed) {
  init_seirconn(m, days, seed)
}

#' @rdname ModelSEIRCONN
#' @export
print.epiworld_seirconn <- function(x, ...) {
  print_seirconn(x)
}

#' #' @rdname ModelSEIRCONN
#' #' @export
#' agents_smallworld.epiworld_seirconn <- function(m, n, k, d, p) {
#'   agents_smallworld_sir(m, n, k, d, p)
#' }

#' @rdname ModelSEIRCONN
#' @export
run.epiworld_seirconn <- function(m) {
  run_seirconn(m)
}
