# result <- suppressWarnings(suppressMessages(invisible({
  
  library(epiworldR)
  
#   library(data.table)
#   library(EpiEstim)
#   library(dplyr)
#   library(ggplot2)
#   library(tidyverse)
#   library(EpiNow2)
# })))


model_seircon=ModelSEIRCONN(name="covid",n=50000,prevalence = 0.001,contact_rate = 20,transmission_rate = 0.5,recovery_rate = 1/7,incubation_days = 3)





# run(model_seircon,ndays=50,seed=1912)
#> _________________________________________________________________________
#> |Running the model...
#> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#> | done.
# summary(model_seircon)
#> ________________________________________________________________________________
#> ________________________________________________________________________________
#> SIMULATION STUDY
#> 
#> Name of the model   : Susceptible-Exposed-Infected-Removed (SEIR) (connected)
#> Population size     : 50000
#> Agents' data        : (none)
#> Number of entities  : 0
#> Days (duration)     : 50 (of 50)
#> Number of viruses   : 1
#> Last run elapsed t  : 473.00ms
#> Last run speed      : 5.28 million agents x day / second
#> Rewiring            : off
#> 
#> Global events:
#>  (none)
#> 
#> Virus(es):
#>  - covid (baseline prevalence: 0.10%)
#> 
#> Tool(s):
#>  (none)
#> 
#> Model parameters:
#>  - Avg. Incubation days : 3.0000
#>  - Contact rate         : 20.0000
#>  - Prob. Recovery       : 0.1429
#>  - Prob. Transmission   : 0.5000
#> 
#> Distribution of the population at time 50:
#>   - (0) Susceptible : 49950 -> 0
#>   - (1) Exposed     :    50 -> 0
#>   - (2) Infected    :     0 -> 140
#>   - (3) Recovered   :     0 -> 49860
#> 
#> Transition Probabilities:
#>  - Susceptible  0.60  0.40  0.00  0.00
#>  - Exposed      0.00  0.66  0.34  0.00
#>  - Infected     0.00  0.00  0.86  0.14
#>  - Recovered    0.00  0.00  0.00  1.00



# plot(model_seircon)



saver <- make_saver(
  "total_hist",
  "transmission",
  "transition",
  "reproductive",
  "generation"
)
run_multiple(model_seircon,ndays=50,nsim=100,seed=1972,saver=saver)
#> Starting multiple runs (100) using 1 thread(s)
#> _________________________________________________________________________
#> _________________________________________________________________________
#> ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.
#>  done.


res1 <- run_multiple_get_results(model_seircon)
#> Warning in run_multiple_get_results(model_seircon): When retrieving the saved
#> results, for the case of transmission, there were no observations.
#> Error in attributes(.Data) <- c(attributes(.Data), attrib): all attributes must have names [4 does not]
res1$reproductive
#> Error in eval(expr, envir, enclos): object 'res1' not found