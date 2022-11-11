 #------------------------------------------------------------------------------
# BUILDING AND INITIALIZING SEIR MODEL
# library(epiworldR)

# sir <- ModelSIR(name = "COVID-19", prevalence = 0.01, infectiousness = 0.9, 
#                 recovery = 0.1)


# # Adding a Small world population
# agents_smallworld(
#   sir,
#   n = 1000,
#   k = 5,
#   d = FALSE,
#   p = .01
#   )

# # Initializing
# init(sir, days = 100, seed = 1912)
# # Running and printing
# # queuing_off(sir)
# run(sir)
# sir


find_scale <- function(x) {
  res <- 10^(floor(log10(x)) + 1 - 3)
  if (res < 1000)
    return(1)
  res
}

  
plot_epi <- function(x, main = "", counts_scale, ...) { 
  x <- get_hist_total(x)
  status_names <- sort(unique(x$status))
  
  # If the user didn't say what scale
  if (missing(counts_scale))
    counts_scale <- find_scale(max(x$counts))
  
  x$counts <- x$counts/counts_scale
  
  # Initialize date vector of size length for status names
  date_candidates <- integer(length = length(status_names)) 
  # Identify max date when the counts stop significantly changing by status
  
  benchmark_value <- diff(range(x$counts))/200 # 0.5% of range
  
  for (i in 1:length(status_names)) {
    date_candidates[i] <- with(
      x[x$status == status_names[i],], 
      sum(abs(diff(counts)) >benchmark_value )
      )
  }
  # Round the maximum date up to the nearest 10th 
  max_date <- ceiling(max(date_candidates) / 10) * 10
  
  # Defining range of x values by max date as the max
  x <- x[x$dates < max_date,]
  # Defining range of y values  
  counts_range <- range(x$counts)

  # Plot the first status
  with(
    x[x$status == status_names[1],], 
    plot(
      x = dates, y = counts, 
      type = 'l', col = 1, ylim = counts_range, 
      xlab = "Days", 
      ylab = ifelse(
        counts_scale == 1,
        "Population",
        paste("Population (", counts_scale, "'s)", sep = "")
        ), 
      main = main
    )
  )
  
  # Plot the remaining statuses
  for (i in 2:length(status_names)) {
    with(
      x[x$status == status_names[i],],
      lines(x = dates, y = counts, type = 'l', col = i)
    )
  }
  # Legend
  legend("right", legend = status_names, col = 1:length(status_names), lty = 1, 
         lwd = 2, bty = "n")
}

plot_epi(sir, main = "SIR Model")
 
