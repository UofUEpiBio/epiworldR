# Test just this file: tinytest::run_test_file("inst/tinytest/test-agents-and-dist-virus-tool.R")

###############################################################################
# Test Group 1: ModelMeaslesQuarantine
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
  distfun = distribute_virus_to_set(10L)
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

###############################################################################
# Fourth case: A random set of two agents are vaccinated
# - At least one susceptible always at the end.
# - If bad luck, the single case doesn't spread (199 susceptible)
# - Most of the time two susceptible.
set_distribution_tool(
  get_tool(abm, 0),
  distfun = distribute_tool_randomly(2L, FALSE)
)

set_distribution_virus(
  get_virus(abm, 0),
  distfun = distribute_virus_to_set(2L)
)

run_multiple(
  abm,
  ndays = 50,
  nsims = 100,
  seed = 1122,
  saver = make_saver("total_hist"),
  nthreads = 2
)

ans <- run_multiple_get_results(abm, nthreads = 2)$total_hist |>
  data.table::as.data.table()

ans <- ans[date == 50 & state == "Susceptible"]
expect_true(ans[, all(counts %in% c(1, 2, 199))])


###############################################################################
# Test Group 2: ModelSEIRMixing
###############################################################################
# First case: Everyone is vaccinated, so only the 
# index case is infected
e1 <- entity("Population 1", 300, FALSE)
e2 <- entity("Population 2", 300, FALSE)
e3 <- entity("Population 3", 300, FALSE)

contact_matrix <- c(
  c(1/3, 1/3, 1/3),
  c(1/3, 1/3, 1/3),
  c(1/3, 1/3, 1/3)
) |> as.double() |> matrix(byrow = TRUE, nrow = 3)

abm <- ModelSEIRMixing(
  name = "Measles",
  n = 900,
  prevalence = 1 / 900,
  contact_rate = 15 / 0.9 / 4,
  transmission_rate = 0.9,
  recovery_rate = 1 / 10,
  incubation_days = 1/12,
  contact_matrix = contact_matrix
)

measles_vaccine <- tool(
  name = "Measles Vaccine",
  susceptibility_reduction = 1.0,
  transmission_reduction = 1.0,
  recovery_enhancer = 1.0,
  death_reduction = 1.0,
  prevalence = 1.0,
  as_proportion = TRUE
)

set_distribution_tool(
  tool = measles_vaccine,
  distfun = distribute_tool_randomly(2L, TRUE)
)

add_tool(abm, measles_vaccine)

abm |>
  add_entity(e1) |>
  add_entity(e2) |>
  add_entity(e3)

# Changing the distribution of the virus
set_distribution_virus(
  get_virus(abm, 0),
  distfun = distribute_virus_to_set(10L)
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
    n = c("5", "895"), state = c("Exposed", "Susceptible")), class = "table")
)

###############################################################################
# Third case: Only agents 25 and 36 are vaccinated:

set_distribution_tool(
  get_tool(abm, 0),
  distfun = distribute_tool_to_set(c(25L, 36L))
)

# Making sure that the agents become infected
set_param(abm, "Prob. Transmission", 1.0)
set_param(abm, "Contact rate", 100)

run(abm, 100, 200)
states <- get_agents_states(abm)

expect_true(
  all(states[c(26, 37)] == "Susceptible") &&
    all(states[-c(26, 37)] == "Recovered")
)

expect_stdout(print(get_agents(abm)), "more agents ...")

###############################################################################
# Fourth case: A random set of two agents are vaccinated
# - At least one susceptible always at the end.
# - If bad luck, the single case doesn't spread (199 susceptible)
# - Most of the time two susceptible.
set_distribution_tool(
  get_tool(abm, 0),
  distfun = distribute_tool_randomly(2L, FALSE)
)

set_distribution_virus(
  get_virus(abm, 0),
  distfun = distribute_virus_to_set(2L)
)

run_multiple(
  abm,
  ndays = 50,
  nsims = 100,
  seed = 1122,
  saver = make_saver("total_hist"),
  nthreads = 2
)

ans <- run_multiple_get_results(abm, nthreads = 2)$total_hist |>
  data.table::as.data.table()

ans <- ans[date == 50 & state == "Susceptible"]
expect_true(ans[, all(counts %in% c(1, 2, 199))])
