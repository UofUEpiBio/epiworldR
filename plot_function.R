#-------------------------------------------------------------------------------
# BUILDING AND INITIALIZING SEIR MODEL
library(epiworldR)
seir <- ModelSEIR(name = "COVID-19", prevalence = 0.01, infectiousness = 0.9, incubation_days = 4, recovery = 0.3)

# Adding a Small world population 
agents_smallworld(
  seir,
  n = 200000,
  k = 5,
  d = FALSE,
  p = .01
  )

# Initializing 
init(seir, days = 100, seed = 1912)

# Running and printing
#queue_off(seir)
run(seir)
seir

#-------------------------------------------------------------------------------
# PLOTTING FUNCTION

#1 Specify dates, status, and counts
test <- get_hist_total(seir)

#2 Build seir plotting function
plot_seir <- function(x){
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

# Run seir plotting function
plot_seir(test)
