library(epiworldR)
library(epiworldRfaster)

library(microbenchmark)


sir <- epiworldR::ModelSEIRCONN(
  name                = "COVID-19",
  prevalence          = 0.01, 
  n                   = 100000,
  contact_rate        = 4, 
  incubation_days     = 7, 
  transmission_rate   = 0.6,
  recovery_rate       = 0.5
) |> 
  epiworldR::add_virus(
    epiworldR::virus("COVID-19-beta", 0.01, 0.6, 0.5, 7), .2
    ) |>
  epiworldR::verbose_off()


sirfast <- epiworldRfaster::ModelSEIRCONN(
  name                = "COVID-19",
  prevalence          = 0.01, 
  n                   = 100000,
  contact_rate        = 4, 
  incubation_days     = 7, 
  transmission_rate   = 0.6,
  recovery_rate       = 0.5
) |> 
  epiworldRfaster::add_virus(
    epiworldRfaster::virus("COVID-19-beta", 0.01, 0.6, 0.5, 7), .5
    ) |>
  epiworldRfaster::verbose_off()

  

res <- microbenchmark(
  old = epiworldR::run(sir, ndays = 100, seed = 1912),
  new = epiworldRfaster::run(sirfast, ndays = 100, seed = 1912),
  times = 10
)

res
summary(sirfast)
summary(sir)

boxplot(res)

