
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R-CMD-check](https://github.com/UofUEpi/epiworldR/actions/workflows/r.yml/badge.svg)](https://github.com/UofUEpi/epiworldR/actions/workflows/r.yml)

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

# Initializing 
init(sir, days = 100, seed = 1912)

# Running and printing
queuing_off(sir)
run(sir)
#> Running the model...
#> _________________________________________________________________________
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
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
#> Last run elapsed t  : 458.00ms
#> Last run speed      : 21.79 million agents x day / second
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
#>  - (0) Susceptible :  99000 -> 2130
#>  - (1) Infected    :   1000 -> 330
#>  - (2) Recovered   :      0 -> 97540
#> 
#> Transition Probabilities:
#>  - Susceptible  0.96  0.04  0.00
#>  - Infected     0.00  0.70  0.30
#>  - Recovered    0.00  0.00  1.00
```

Visualizing the outputs

``` r
plot(sir)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
