 # #-------------------------------------------------------------------------------
 #  # BUILDING AND INITIALIZING SEIR MODEL
 #  library(epiworldR)
 #  seir <- ModelSEIR(name = "COVID-19", prevalence = 0.01, infectiousness = 0.9, incubation_days = 4, recovery = 0.3)
 #  
 #  # Adding a Small world population 
 #  agents_smallworld(
 #    seir,
 #    n = 200000,
 #    k = 5,
 #    d = FALSE,
 #    p = .01
 #    )
 #  
 #  # Initializing 
 #  init(seir, days = 100, seed = 1912)
 #  
 #  # Running and printing
 #  queuing_off(seir)
 #  run(seir)
 #  seir

plot_epi <- function(x, main = "", ...) { 
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
                                       main = main))
  # Plot the remaining statuses
  for (i in 2:length(vnames)) {
    with(x[x$status == vnames[i],] ,lines(x = dates, y = counts, type = 'l', 
                                          col = i))
  }
  # Legend
  legend("right", legend = vnames, col = 1:length(vnames), lty = 1, lwd = 2, 
         bty = "n")
}

# plot_epi(seir, main = "SEIR Model")


