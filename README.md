
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
  name           = "COVID-19",
  prevalence     = .01,
  infectiousness = .9,
  recovery       = .3
  ) |>
  # Adding a Small world population 
  agents_smallworld(n = 100000, k = 10, d = FALSE, p = .01) |>
  # Running the model for 100 days
  run(ndays = 100, seed = 1912)
#> _________________________________________________________________________
#> Running the model...
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#>  done.

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
#> Last run elapsed t  : 308.00ms
#> Last run speed      : 32.46 million agents x day / second
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
#>  - Infected     0.00  0.69  0.31
#>  - Recovered    0.00  0.00  1.00
```

Visualizing the outputs

``` r
plot(sir)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
