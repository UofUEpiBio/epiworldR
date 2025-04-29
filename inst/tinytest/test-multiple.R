# result <- suppressWarnings(suppressMessages(invisible({
  
n <- 5e3
days <- 50
nsims <- 50

model_seircon=ModelSEIRCONN(
  name="covid",
  n = n,
  prevalence = 0.001,
  contact_rate = 20,
  transmission_rate = 0.5,
  recovery_rate = 1/7,
  incubation_days = 3
  )

run(model_seircon,ndays=days,seed=1912)

saver <- make_saver(
  "total_hist",
  "transmission",
  "transition",
  "reproductive",
  "generation"
)

run_multiple(
  model_seircon,
  ndays=days,
  nsims=nsims,
  seed=1972,
  saver=saver,
  nthreads = 2L
  )

res1 <- run_multiple_get_results(model_seircon, nthreads = 2L)

res1 <- lapply(res1, data.table::as.data.table)

avg_infected <- res1$total_hist[date==50 & state=="Infected"]$counts |> mean()
avg_recovered <- res1$total_hist[date==50 & state=="Recovered"]$counts |> mean()

expect_true(abs(avg_infected - 13.8) < 50)
expect_true(abs(avg_recovered - 4986) < 50)
