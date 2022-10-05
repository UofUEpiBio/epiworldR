
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epiworldR

<!-- badges: start -->
<!-- badges: end -->

This R package is a wrapper of the C++ library
[epiworld](https://github.com/UofUEpi/epiworld). It provides a general
framework for modeling disease transmission using Agent-Based Models.

## Installation

You can install the development version of epiworldR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("UofUEpi/epiworldR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(epiworldR)

# Creating a SIR model
sir <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  infectiousness = .9,
  recovery = .3
  )

# Adding a Small world population 
agents_smallworld(
  sir,
  n = 100000,
  k = 5,
  d = FALSE,
  p = .01
  )
#> [1] 0

# Initializing 
init(sir, days = 100, seed = 1912)
#> [1] 0

# Running and printing
run(sir)
#> Running the model...
#> _________________________________________________________________________
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#> [1] 0
sir
#> 
#> ________________________________________________________________________________
#> SIMULATION STUDY
#> 
#> Name of the model   : Susceptible-Infected-Recovered (SIR)
#> Population size     : 100000
#> Number of entitites : 0
#> Days (duration)     : 100 (of 100)
#> Number of variants  : 1
#> Last run elapsed t  : 1.00s
#> Last run speed      : 6.99 million agents x day / second
#> Rewiring            : off
#> 
#> Virus(es):
#>  - COVID-19 (baseline prevalence: 1.00%)
#> 
#> Tool(s):
#>  (none)
#> 
#> Model parameters:
#>  - Infectiousness    : 0.9000
#>  - Prob. of Recovery : 0.3000
#> 
#> Distribution of the population at time 100:
#>  - (0) Susceptible :  99000 -> 0
#>  - (1) Infected    :   1000 -> 0
#>  - (2) Recovered   :      0 -> 100000
#> 
#> Transition Probabilities:
#>  - Susceptible  0.87  0.13  0.00
#>  - Infected     0.00  0.70  0.30
#>  - Recovered    0.00  0.00  1.00
```

Visualizing the outputs

``` r
x <- get_hist_total(sir)
x$counts <- x$counts/1000
x <- x[x$dates < 50,]

with(
  x[x$status == "Susceptible",],
  plot(
    x = dates, y = counts, type = "l", col = "blue", ylim = range(x$counts),
    ylab = "Population (thousands)", xlab = "days", main = "SIR model")
  )

with(
  x[x$status == "Infected",],
  lines(x = dates, y = counts, col = "red")
  )

with(
  x[x$status == "Recovered",],
  lines(x = dates, y = counts, col = "darkgreen")
  )

legend(
  "right",
  legend = c("Susceptible", "Infected", "Removed"),
  col    = c("blue", "red", "darkgreen"),
  lty    = 1,
  lwd    = 2,
  bty    = "n"
  )
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
