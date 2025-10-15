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
  ) |>
  add_tool(
    tool = tool(
      name = "test_tool",
      prevalence = .5,
      as_proportion = TRUE,
      susceptibility_reduction = 0.0,
      transmission_reduction = 0.0,
      recovery_enhancer = 0.0,
      death_reduction = 0.0
      )
  )

run(model_seircon,ndays=days,seed=1912)

saver <- make_saver(
  "total_hist",
  "virus_hist",
  "virus_info",
  "tool_info",
  "tool_hist",
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

res1 <- run_multiple_get_results(
  model_seircon, nthreads = 1L, freader = data.table::fread
  )


