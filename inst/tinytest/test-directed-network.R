# Directed networks -------------------------------------------------------------
# A directed tie i -> j is kept by its source: j is one of i's neighbors, and
# agents catch the virus from their own neighbors, so j can infect i but not
# the reverse. Along the path 0 -> 1 -> ... -> 9, a virus seeded at 9 reaches
# everyone and one seeded at 0 reaches no one.
n <- 10L
model <- ModelSIR("a virus", prevalence = 0, transmission_rate = 1,
                  recovery_rate = 0)
agents_from_edgelist(
  model, source = 0L:(n - 2L), target = 1L:(n - 1L), size = n, directed = TRUE
)
verbose_off(model)

set_distribution_virus(get_virus(model, 0), distribute_virus_to_set(n - 1L))
run(model, ndays = 20, seed = 1)
expect_equal(sum(get_agents_states(model) == "Infected"), n)

set_distribution_virus(get_virus(model, 0), distribute_virus_to_set(0L))
run(model, ndays = 20, seed = 1)
expect_equal(sum(get_agents_states(model) == "Infected"), 1L)

# The ties are kept as given, not mirrored
net <- get_network(model)
expect_equal(nrow(net), n - 1L)
expect_equal(net$from, 0L:(n - 2L))
expect_equal(net$to, 1L:(n - 1L))
