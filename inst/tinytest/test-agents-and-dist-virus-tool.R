
###############################################################################
# First case: Everyone is vaccinated, so only the 
# index case is infected
abm <- ModelMeaslesQuarantine(
  n = 200,
  prop_vaccinated = 1.0,
  vax_efficacy = 1.0
)

# Changing the distribution of the virus
set_distribution_virus(
  get_virus(abm, 0),
  distfun = distribute_virus_set(10L)
)

run(abm, 100, 200)

states <- get_agents_states(abm)
expect_true(
  (states[11] != "Susceptible") && all(states[-11] == "Susceptible")
)

###############################################################################
# Second case: Distribute the virus randomly to 5 agents
set_distribution_virus(
  get_virus(abm, 0),
  distfun = distribute_virus_randomly(5L, FALSE)
)

run_multiple(
  abm,
  ndays = 50,
  nsims = 100,
  seed = 1122,
  saver = make_saver("total_hist"),
  nthreads = 2
)

ans <- run_multiple_get_results(abm, nthreads = 2)
ans <- data.table::as.data.table(ans$total_hist)[date == 0]
ans <- ans[, .(n = sum(counts)), by = .(sim_num, state)][n > 0]
expect_equal(
  ans[, table(n, state)],
  structure(c(100L, 0L, 0L, 100L), dim = c(2L, 2L), dimnames = list(
    n = c("5", "195"), state = c("Exposed", "Susceptible")), class = "table")
)


###############################################################################
# Third case: Only agents 24 and 36 are vaccinated:

set_distribution_tool(
  get_tool(abm, 0),
  distfun = distribute_tool_to_set(c(25L, 36L))
)

# Making sure that the agents become infected
set_param(abm, "Transmission rate", 1.0)
set_param(abm, "Contact rate", 100)

run(abm, 100, 200)
states <- get_agents_states(abm)

expect_true(
  all(states[c(26, 37)] == "Susceptible") &&
  all(states[-c(26, 37)] == "Recovered")
)

expect_stdout(print(get_agents(abm)), "more agents ...")
