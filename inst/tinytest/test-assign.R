# Placement: viruses, tools, initial stages and groups (assign_*()) ----------
# Agent IDs are 0-based.

vaccine <- tool(
  "Vaccine",
  susceptibility_reduction = 1, transmission_reduction = 0,
  recovery_enhancer = 0, death_reduction = 0
)

# A subset to draw from is a group of a mixing layer: 30 of agents 900:999
# get the virus, and half of them the tool
abm <- ModelSIR(
  "a virus", n = 1000, prevalence = 0,
  transmission_prob = 0, recovery_rate = 0
)
add_layer(
  abm,
  layer_mixing(
    "subset",
    groups       = c(rep(NA, 900), rep(0, 100)),
    contact_prob = 0
  )
)
assign_virus(abm, 0, assign_group_counts("subset", 30))
add_tool(abm, vaccine, assign_group_shares("subset", .5))
run(abm, ndays = 0, seed = 123)

infected <- which(get_agents_stages(abm) == "Infected") - 1L
expect_equal(length(infected), 30L)
expect_true(all(infected %in% 900:999))
expect_equal(sum(get_hist_tool(abm)$counts), 50)

# The draws change with the seed, and repeat with it
run(abm, ndays = 0, seed = 123)
expect_identical(which(get_agents_stages(abm) == "Infected") - 1L, infected)
run(abm, ndays = 0, seed = 321)
expect_false(identical(which(get_agents_stages(abm) == "Infected") - 1L, infected))

# The seeds start in the infection's target (a carrier stage), and
# set_initial() places agents without the virus
abm <- ModelSIR(
  "a virus", n = 1000, prevalence = 0,
  transmission_prob = 0, recovery_rate = 0
)
add_param(abm, "Contact rate", 0)
add_layer(
  abm,
  layer_mixing("everyone", groups = rep(0, 1000), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)
assign_virus(abm, 0, assign_count(100))
set_initial(abm, assign_counts(c(Recovered = 300)))
run(abm, ndays = 0, seed = 1)
expect_equal(
  as.vector(table(factor(get_agents_stages(abm), get_stages(abm)))),
  c(600, 100, 300)
)

# A stage that cannot hold the virus is not a seed's stage, and labels are
# checked against the stages
expect_error(
  assign_virus(abm, 0, assign_counts(c(Recovered = 10))) |> run(0, 1),
  "carrier|virus"
)
expect_error(
  set_initial(abm, assign_shares(c(Recover = .1))),
  "no stage \"Recover\""
)

# Fixed labels, and labels from a function (once per agent and run)
assign_virus(abm, 0, assign_labels(c(rep("Infected", 5), rep(NA, 995))))
set_initial(abm, assign_by(function(agent, model) NA))
run(abm, ndays = 0, seed = 1)
expect_equal(which(get_agents_stages(abm) == "Infected") - 1L, 0:4)
expect_true(isTRUE(attr(abm, "single_threaded")))

# Invalid placements
expect_error(assign_share(1.5), "\\[0, 1\\]")
expect_error(assign_count(2.5), "whole")
expect_error(assign_agents(-1), "agent ids")
expect_error(add_tool(abm, vaccine, .5), "assign_")

# assign_agents() with a large set ---------------------------------------------
# Placing k agents takes about as long as a random share (it once took O(k^2)).
n <- 100000L
ids <- seq(0L, n - 1L, by = 2L)

abm <- ModelSIR(
  "a virus", n = n, prevalence = 0,
  transmission_prob = 0, recovery_rate = 0
)
add_param(abm, "Contact rate", 0)
add_layer(
  abm,
  layer_mixing("everyone", groups = rep(0, n), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)
assign_virus(abm, 0, assign_agents(ids))
add_tool(abm, vaccine, assign_agents(ids))

elapsed <- system.time(run(abm, ndays = 0, seed = 123))[["elapsed"]]

expect_equal(which(get_agents_stages(abm) == "Infected") - 1L, ids)
expect_equal(sum(get_hist_tool(abm)$counts), length(ids))
expect_true(elapsed < 2)
