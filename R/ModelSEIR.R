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
plot.epiworld_seir <- function(x, ...) { # col = NULL
  x <- get_hist_total(x)
  vnames <- sort(unique(x$status)) 
  x$counts <- x$counts/1000
  x <- x[x$dates < 50,]
  counts_range <- range(x$counts)

  # Plot the first status
  with(x[x$status == vnames[1],], plot(x = dates, y = counts, 
                                       type = 'l', col = 1, ylim = counts_range, 
                                       xlab = "Days", 
                                       ylab = "Population (thousands)", 
                                       main = "SEIR Model"))

  # Plot the remaining statuses
  for (i in 2:length(vnames)) {
    with(x[x$status == vnames[i],] ,lines(x = dates, y = counts, type = 'l', 
                                          col = i))
  }

  legend("right", legend = vnames, col = 1:length(vnames), lty = 1, lwd = 2, 
         bty = "n")
}
