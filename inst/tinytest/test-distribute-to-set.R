# Tools and viruses distributed to a set of agents -----------------------------
# Agent IDs are 0-based.

vaccine <- tool(
  "Vaccine", prevalence = 0, as_proportion = FALSE,
  susceptibility_reduction = 1, transmission_reduction = 0,
  recovery_enhancer = 0, death_reduction = 0
)

# distribute_virus_randomly() with agents_ids samples within the set
abm <- ModelSIRCONN(
  "a virus", n = 1000, prevalence = 0, contact_rate = 0,
  transmission_rate = 0, recovery_rate = 0
)
set_distribution_virus(
  get_virus(abm, 0),
  distribute_virus_randomly(30, as_proportion = FALSE, agents_ids = 900:999)
)
run(abm, ndays = 0, seed = 123)

infected <- which(get_agents_states(abm) == "Infected") - 1L
expect_equal(length(infected), 30L)
expect_true(all(infected %in% 900:999))

# distribute_tool_randomly() with agents_ids samples within the set (it used
# to pick among agents 0..(n-1), n being the size of the set). Everyone but
# the vaccinated gets infected, so the susceptible are the tool holders.
abm <- ModelSIRCONN(
  "a virus", n = 1000, prevalence = 0, contact_rate = 20,
  transmission_rate = 1, recovery_rate = 0
)
set_distribution_virus(get_virus(abm, 0), distribute_virus_to_set(0:9))
add_tool(abm, vaccine)

for (d in list(
  list(fun = distribute_tool_randomly(40, FALSE, agents_ids = 900:999), n = 40),
  list(fun = distribute_tool_randomly(0.5, TRUE, agents_ids = 900:999), n = 50)
)) {

  set_distribution_tool(get_tool(abm, 0), d$fun)
  run(abm, ndays = 30, seed = 123)

  vaccinated <- which(get_agents_states(abm) == "Susceptible") - 1L
  expect_equal(length(vaccinated), d$n)
  expect_true(all(vaccinated %in% 900:999))

}

# distribute_*_to_set() with a large set --------------------------------------
# The set used to be copied every time an agent received the tool/virus, so
# distributing to k agents took O(k^2): several seconds for k = 50,000. It now
# takes about as long as distribute_*_randomly().
n <- 100000L
ids <- seq(0L, n - 1L, by = 2L)

abm <- ModelSIRCONN(
  "a virus", n = n, prevalence = 0, contact_rate = 0,
  transmission_rate = 0, recovery_rate = 0
)
set_distribution_virus(get_virus(abm, 0), distribute_virus_to_set(ids))
add_tool(abm, vaccine)
set_distribution_tool(get_tool(abm, 0), distribute_tool_to_set(ids))

elapsed <- system.time(run(abm, ndays = 0, seed = 123))[["elapsed"]]

expect_equal(which(get_agents_states(abm) == "Infected") - 1L, ids)
expect_equal(sum(get_hist_tool(abm)$counts), length(ids))
expect_true(elapsed < 2)
