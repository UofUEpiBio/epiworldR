# library(epiworldRdev)
# library(epiworldR)

library(microbenchmark)

ns <- c(1e3, 5e3, 1e4, 5e4, 1e5, 5e5)
ans <- vector("list", length(ns))
names(ans) <- as.character(ns)
for (n in ns) {

  sir <- epiworldR::ModelSEIR(
    name                = "COVID-19",
    prevalence          = 0.01,
    incubation_days     = 7,
    transmission_rate   = 0.6,
    recovery_rate       = 0.5
  ) |>
    epiworldR::agents_smallworld(n = n, k = 20, p = 0.0, d = FALSE) |>
    epiworldR::add_virus(
      epiworldR::virus("COVID-19-beta", 0.01, 0.6, 0.5, 7), .2
    ) |>
    epiworldR::verbose_off()


  sirfast <- epiworldRdev::ModelSEIR(
    name                = "COVID-19",
    prevalence          = 0.01,
    incubation_days     = 7,
    transmission_rate   = 0.6,
    recovery_rate       = 0.5
  ) |>
    epiworldRdev::agents_smallworld(n = n, k = 20, p = 0, d = FALSE) |>
    epiworldRdev::add_virus(
      epiworldRdev::virus("COVID-19-beta", 0.01, 0.6, 0.5, 7), .2
    ) |>
    epiworldRdev::verbose_off()


  ans[[as.character(n)]] <- microbenchmark(
    old = epiworldR::run(sir, ndays = 100, seed = 1912),
    new = epiworldRdev::run(sirfast, ndays = 100, seed = 1912),
    times = 10
  )

  message("Simulation with ", n, " individuals finished.")

}

saveRDS(ans, "playground/benchmark-seir.rds")
