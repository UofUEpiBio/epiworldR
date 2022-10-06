#' Susceptible Infected Susceptible model (SEIR)
#'
#' @param name Name of the virus
#' @param prevalence a number
#' @param infectiousness a number
#' @param incubation_days a number
#' @param recovery a number
#' @param x to be documented
#' @param ... to be documented
#' @export
#' @family Models
#' @aliases epiworld_seir
ModelSEIR <- function(
    name, prevalence, infectiousness, incubation_days, recovery
) {
  
  structure(
    ModelSEIR_cpp(name, prevalence, infectiousness, incubation_days, recovery),
    class = "epiworld_seir"
  )
  
}

#' @param m to be documented
#'
#' @param days to be documented
#' @param seed to be documented
#'
#' @rdname ModelSEIR
#' @export
init.epiworld_seir <- function(m, days, seed) {
  init_sir(m, days, seed)
}

#' @rdname ModelSEIR
#' @export
print.epiworld_seir <- function(x, ...) {
  print_sir(x)
}

#' @param n to be documented
#' @param k to be documented
#' @param d to be documented
#' @param p to be documented
#'
#' @rdname ModelSEIR
#' @export
agents_smallworld.epiworld_seir <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSEIR
#' @export
run.epiworld_seir <- function(m) {
  run_sir(m)
}

#' @rdname ModelSEIR
#' @export
plot.epiworld_seir <- function(seir_model) {
  x <- get_hist_total(seir_model)
    x$counts <- x$counts/1000
    x <- x[x$dates < 50,]
    
    with(
      x[x$status == "Susceptible",],
      plot(
        x = dates, y = counts, type = "l", col = "blue", ylim = range(x$counts),
        ylab = "Population (thousands)", xlab = "days", main = "SEIR model")
      )
    
    with(
      x[x$status == "Exposed",],
      lines(x = dates, y = counts, col = "purple")
      )
    
    with(
      x[x$status == "Infected",],
      lines(x = dates, y = counts, col = "red")
      )
    
    with(
      x[x$status == "Removed",],
      lines(x = dates, y = counts, col = "darkgreen")
      )
    
    legend(
      "right",
      legend = c("Susceptible", "Exposed", "Infected", "Removed"),
      col    = c("blue", "purple", "red", "darkgreen"),
      lty    = 1,
      lwd    = 2,
      bty    = "n"
      )
}
