# result <- suppressWarnings(suppressMessages(invisible({
  
n <- 5e3
days <- 50
nsims <- 50

model_seircon <- ModelSEIR(
  name="covid",
  n = n,
  prevalence = 0.001,
  transmission_prob = 0.5,
  recovery_rate = 1/7,
  latent_period = 3
  )
add_param(model_seircon, "Contact rate", 20)
add_layer(
  model_seircon,
  layer_mixing("everyone", groups = rep(0, n), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
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
avg_recovered <- res1$total_hist[date==50 & state=="Removed"]$counts |> mean()

expect_true(abs(avg_infected - 13.8) < 50)
expect_true(abs(avg_recovered - 4986) < 50)

# run_multiple() should enforce single-thread execution when globalevent_fun() is used
model_single_thread <- ModelSIR(
  name = "covid",
  n = 300,
  prevalence = 0.02,
  transmission_prob = 0.4,
  recovery_rate = 0.2
)
add_param(model_single_thread, "Contact rate", 4)
add_layer(
  model_single_thread,
  layer_mixing("everyone", groups = rep(0, 300), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

verbose_off(model_single_thread)
add_globalevent(
  model_single_thread,
  globalevent_fun(function(m) invisible(NULL), "No-op event")
)

expect_true(isTRUE(attr(model_single_thread, "single_threaded")))

saver_single_thread <- make_saver("total_hist")
expect_warning(
  run_multiple(
    model_single_thread,
    ndays = 10,
    nsims = 4,
    seed = 42,
    saver = saver_single_thread,
    verbose = FALSE,
    nthreads = 2L
  ),
  "nthreads"
)

results_single_thread <- run_multiple_get_results(
  model_single_thread,
  nthreads = 1L
)
expect_true(nrow(results_single_thread$total_hist) > 0)
