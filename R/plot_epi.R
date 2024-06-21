 #------------------------------------------------------------------------------
# BUILDING AND INITIALIZING SEIR MODEL
# library(epiworldR)

# sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, transmission_rate = 0.9, 
#                 recovery_rate = 0.1)


# # Adding a Small world population
# agents_smallworld(
#   sir,
#   n = 1000,
#   k = 5,
#   d = FALSE,
#   p = .01
#   )

# # Running and printing
# # queuing_off(sir)
# run(sir, ndays = 100, seed = 1912)
# sir


find_scale <- function(x) {
  res <- 10^(floor(log10(x)) + 1 - 3)
  if (res < 1000)
    return(1)
  res
}

#' @noRd
#' @importFrom graphics legend
plot_epi <- function(x, main = "", counts_scale, ...) UseMethod("plot_epi")

#' @export
plot_epi.epiworld_model <- function(
    x, main = "",
    counts_scale,
    ...
) {
  
  plot_epi(
    x = get_hist_total(x),
    main = main,
    counts_scale = counts_scale,
    ...
  )
  
}

#' @export
plot_epi.epiworld_hist_virus <- function(
    x, main = "",
    counts_scale,
    ...
) {
  
  res <- lapply(sort(unique(x$id)), function(i) x[x$id == i,])
  
  lapply(res, function(r) {
    plot_epi.epiworld_hist(r, main = paste0("Variant id ", r$id[1]))
    })
  invisible(x)
  
}

#' @export
plot_epi.epiworld_hist <- function(
    x, main = "",
    counts_scale,
    ...
    ) { 
  
  curves      <- x
  state_names <- attr(curves, "states")
  
  # If the user didn't say what scale
  if (missing(counts_scale))
    counts_scale <- find_scale(max(curves$counts))
  
  curves$counts <- curves$counts/counts_scale
  
  # Initialize date vector of size length for state names
  date_candidates <- integer(length = length(state_names)) 
  # Identify max date when the counts stop significantly changing by state
  
  benchmark_value <- diff(range(curves$counts))/200 # 0.5% of range
  
  for (i in 1L:length(state_names)) {
    date_candidates[i] <- with(
      curves[curves$state == state_names[i],], 
      sum(abs(diff(counts)) > benchmark_value )
      )
  }
  # Round the maximum date up to the nearest 10th 
  max_date <- min(
    diff(range(curves$date)),
    max(ceiling(max(date_candidates) / 10L) * 10L, 10L)
  )
  
  # Defining range of x values by max date as the max
  curves <- curves[curves$date < max_date,]
  # Defining range of y values  
  counts_range <- range(curves$counts)

  # Plot the first state
  with(
    curves[curves$state == state_names[1L],], 
    graphics::plot(
      x    = date,
      y    = counts, 
      type = 'l', 
      col  = 1,
      ylim = counts_range, 
      xlab = "Day (step)", 
      ylab = ifelse(
        counts_scale == 1L,
        "Population",
        paste("Population (", counts_scale, "'s)", sep = "")
        ), 
      main = main
    )
  )
  
  # Plot the remaining states
  for (i in 2L:length(state_names)) {
    
    with(
      curves[curves$state == state_names[i],],
      graphics::lines(
        x    = date,
        y    = counts,
        type = 'l',
        col  = i
        )
    )
    
  }
  
  # Legend
  graphics::legend(
    "right",
    legend = state_names,
    col    = 1L:length(state_names),
    lty    = 1L, 
    lwd    = 2L,
    bty    = "n"
    )
}

# plot_epi(sir, main = "SIR Model")
 
