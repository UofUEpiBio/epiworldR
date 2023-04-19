
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

# Examples

## SIR model using a random graph

This is a basic example which shows you how to solve a common problem:

``` r
library(epiworldR)

# Creating a SIR model
sir <- ModelSIR(
  name           = "COVID-19",
  prevalence     = .01,
  infectiousness = .7,
  recovery       = .3
  ) |>
  # Adding a Small world population 
  agents_smallworld(n = 100000, k = 10, d = FALSE, p = .01) |>
  # Running the model for 50 days
  run(ndays = 50, seed = 1912)
#> _________________________________________________________________________
#> |Running the model...
#> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#> | done.

sir
#> 
#> ________________________________________________________________________________
#> SIMULATION STUDY
#> 
#> Name of the model   : Susceptible-Infected-Recovered (SIR)
#> Population size     : 100000
#> Number of entities  : 0
#> Days (duration)     : 50 (of 50)
#> Number of variants  : 1
#> Last run elapsed t  : 177.00ms
#> Last run speed      : 28.18 million agents x day / second
#> Rewiring            : off
#> 
#> Virus(es):
#>  - COVID-19 (baseline prevalence: 1.00%)
#> 
#> Tool(s):
#>  (none)
#> 
#> Model parameters:
#>  - Infectiousness    : 0.7000
#>  - Prob. of Recovery : 0.3000
#> 
#> Distribution of the population at time 50:
#>  - (0) Susceptible :  99000 -> 822
#>  - (1) Infected    :   1000 -> 415
#>  - (2) Recovered   :      0 -> 98763
#> 
#> Transition Probabilities:
#>  - Susceptible  0.91  0.09  0.00
#>  - Infected     0.00  0.70  0.30
#>  - Recovered    0.00  0.00  1.00
```

Visualizing the outputs

``` r
plot(sir)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## SEIR Model with a fully connected graph

``` r
model_seirconn <- ModelSEIRCONN(
  name                = "COVID-19",
  prevalence          = 0.01, 
  n                   = 10000,
  contact_rate        = 4, 
  incubation_days     = 7, 
  prob_transmission   = 0.6,
  prob_recovery       = 0.5
)

set.seed(132)
run(model_seirconn, ndays = 100)
#> _________________________________________________________________________
#> Running the model...
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#>  done.
model_seirconn
#> 
#> ________________________________________________________________________________
#> SIMULATION STUDY
#> 
#> Name of the model   : Susceptible-Exposed-Infected-Removed (SEIR) (connected)
#> Population size     : 10000
#> Number of entities  : 0
#> Days (duration)     : 100 (of 100)
#> Number of variants  : 1
#> Last run elapsed t  : 26.00ms
#> Last run speed      : 38.27 million agents x day / second
#> Rewiring            : off
#> 
#> Virus(es):
#>  - COVID-19 (baseline prevalence: 1.00%)
#> 
#> Tool(s):
#>  (none)
#> 
#> Model parameters:
#>  - Avg. Incubation days : 7.0000
#>  - Contact rate         : 4.0000
#>  - Prob. Recovery       : 0.5000
#>  - Prob. Transmission   : 0.6000
#> 
#> Distribution of the population at time 100:
#>  - (0) Susceptible :  9900 -> 98
#>  - (1) Exposed     :   100 -> 0
#>  - (2) Infected    :     0 -> 0
#>  - (3) Recovered   :     0 -> 9902
#> 
#> Transition Probabilities:
#>  - Susceptible  0.96  0.04  0.00  0.00
#>  - Exposed      0.00  0.87  0.13  0.00
#>  - Infected     0.00  0.00  0.49  0.51
#>  - Recovered    0.00  0.00  0.00  1.00
```

Computing some key statistics

``` r
plot(model_seirconn)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r

repnum <- get_reproductive_number(model_seirconn)
plot(repnum, type = "b")
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />
