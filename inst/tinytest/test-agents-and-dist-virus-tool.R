# Test just this file: tinytest::run_test_file("inst/tinytest/test-agents-and-dist-virus-tool.R")

###############################################################################
# Test Group 1: ModelSEIR with random mixing
###############################################################################
# First case: Everyone is vaccinated, so only the 
# index case is infected
abm <- ModelSEIR(
  name = "COVID-19",
  n = 200,
  prevalence = 0,
  transmission_prob = 0.9,
  latent_period = 7,
  recovery_rate = 0.3
)
add_param(abm, "Contact rate", 2)
add_layer(
  abm,
  layer_mixing("everyone", groups = rep(0, 200), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)

# Add a vaccine tool that provides 100% protection
vaccine <- tool(
  name = "Vaccine",
  susceptibility_reduction = 1.0,
  transmission_reduction = 1.0,
  recovery_enhancer = 1.0,
  death_reduction = 1.0
)

add_tool(abm, vaccine, assign_share(1))

# Changing the placement of the virus
assign_virus(abm, 0, assign_agents(10L))

run(abm, 100, 200)

states <- get_agents_stages(abm)
expect_true(
  (states[11] != "Susceptible") && all(states[-11] == "Susceptible")
)

###############################################################################
# Second case: Distribute the virus randomly to 5 agents
assign_virus(abm, 0, assign_count(5))

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
ans <- ans[, .(n = sum(counts)), by = .(sim_num, stage)][n > 0]
# All simulations should have 5 exposed and 195 susceptible at start
expect_equal(
  ans[, table(n, stage)],
  structure(c(100L, 0L, 0L, 100L), dim = c(2L, 2L), dimnames = list(
    n = c("5", "195"), stage = c("Exposed", "Susceptible")), class = "table")
)


###############################################################################
# Third case: Only agents 25 and 36 are vaccinated:

assign_tool(abm, 0, assign_agents(c(25L, 36L)))

# Making sure that the agents become infected
set_param(abm, "Transmission probability", 1.0)
set_param(abm, "Contact rate", 100)

run(abm, 100, 200)
states <- get_agents_stages(abm)

# Agents 26 and 37 (indices 25 and 36) should stay susceptible (vaccinated)
# All others should be removed
expect_true(
  all(states[c(26, 37)] == "Susceptible") &&
  all(states[-c(26, 37)] == "Removed")
)

expect_stdout(print(get_agents(abm)), "more agents ...")

###############################################################################
# Fourth case: A random set of two agents are vaccinated
# - At least one susceptible always at the end.
# - If bad luck, the single case doesn't spread (199 susceptible)
# - Most of the time two susceptible.
assign_tool(abm, 0, assign_count(2))

assign_virus(abm, 0, assign_agents(2L))

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

ans <- ans[date == 50 & stage == "Susceptible"]
expect_true(ans[, all(counts %in% c(1, 2, 199))])


###############################################################################
# Test Group 2: ModelSEIR on a mixing layer
###############################################################################
# First case: Everyone is vaccinated, so only the 
# index case is infected
contact_matrix <- (c(
  c(1/3, 1/3, 1/3),
  c(1/3, 1/3, 1/3),
  c(1/3, 1/3, 1/3)
) * (15/0.9/4)) |> as.double() |> matrix(byrow = TRUE, nrow = 3)

abm <- ModelSEIR(
  name = "COVID-19",
  n = 900,
  prevalence = 1 / 900,
  transmission_prob = 0.9,
  recovery_rate = 1 / 10,
  latent_period = 1/12
)

# Everyone mixes evenly
add_layer(abm, layer_mixing("mixing", rep(0, 900), contact_matrix = sum(contact_matrix[1, ])))

vaccine <- tool(
  name = "Vaccine",
  susceptibility_reduction = 1.0,
  transmission_reduction = 1.0,
  recovery_enhancer = 1.0,
  death_reduction = 1.0
)

add_tool(abm, vaccine, assign_share(1))

# Changing the placement of the virus
assign_virus(abm, 0, assign_agents(10L))

run(abm, 100, 200)

states <- get_agents_stages(abm)
expect_true(
  (states[11] != "Susceptible") && all(states[-11] == "Susceptible")
)

###############################################################################
# Second case: Distribute the virus randomly to 5 agents
assign_virus(abm, 0, assign_count(5))

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
ans <- ans[, .(n = sum(counts)), by = .(sim_num, stage)][n > 0]

expect_equal(
  ans[, table(n, stage)],
  structure(c(100L, 0L, 0L, 100L), dim = c(2L, 2L), dimnames = list(
    n = c("5", "895"), stage = c("Exposed", "Susceptible")), class = "table")
)

###############################################################################
# Third case: Only agents 25 and 36 are vaccinated:

assign_tool(abm, 0, assign_agents(c(25L, 36L)))

# Making sure that the agents become infected
set_param(abm, "Transmission probability", 1.0)

run(abm, 100, 200)
states <- get_agents_stages(abm)

expect_true(
  all(states[c(26, 37)] == "Susceptible") &&
    all(states[-c(26, 37)] == "Removed")
)

expect_stdout(print(get_agents(abm)), "more agents ...")

###############################################################################
# Fourth case: A random set of two agents are vaccinated
# - At least one susceptible always at the end.
# - If bad luck, the single case doesn't spread (199 susceptible)
# - Most of the time two susceptible.
assign_tool(abm, 0, assign_count(2))

assign_virus(abm, 0, assign_agents(2L))

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

ans <- ans[date == 50 & stage == "Susceptible"]
expect_true(ans[, all(counts %in% c(1, 2, 199))])

###############################################################################
# Test Group 3: giving agents the model's virus and tool from a global event
###############################################################################
# get_agents(model)[i] is the model's agent (not a copy), so the agent, the
# totals and the tool history agree
abm <- ModelSIR("flu", n = 100, prevalence = 0, transmission_prob = 0,
  recovery_rate = 0)
verbose_off(abm)
add_tool(abm, tool("Vaccine", 0, 0, 0, 0))
give <- function(model) {
  add_virus_agent(get_agents(model)[2], model, get_virus(model, 0))
  add_tool_agent(get_agents(model)[5], model, get_tool(model, 0))
}
add_globalevent(abm, globalevent_fun(give, "give", day = 1))
run(abm, ndays = 3, seed = 1)

stages <- get_agents_stages(abm)
expect_identical(which(stages == "Infected") - 1L, 2L)
expect_equal(
  subset(get_today_total(abm), stage == "Infected")$counts, 1
)
expect_equal(get_stage(abm, get_agents(abm)[2]), "Infected")
expect_equal(get_status(abm, get_agents(abm)[2]), "Free")
expect_equal(
  sum(subset(get_hist_tool(abm), date == 3)$counts), 1
)
