# Contact layers: network layers from generators and edgelists, the layer
# objects shared with the model, mixing layers (across and within groups), and
# the defunct functions.

# A model without layers has no contacts ---------------------------------------
sir <- ModelSIR("flu", n = 1000, prevalence = .05, transmission_prob = 1,
  recovery_rate = .1)
verbose_off(sir)
run(sir, ndays = 20, seed = 1)
expect_equal(nrow(subset(get_transmissions(sir), source >= 0)), 0L)

# A redrawn network: new graph every run, same graph for the same seed ---------
contacts <- layer_network("contacts", graph_smallworld(k = 4, p = .2))
expect_inherits(contacts, "epiworld_layer_network")
add_layer(sir, contacts)
expect_error(add_layer(sir, contacts), "model")

run(sir, ndays = 20, seed = 1)
el1 <- get_edgelist(contacts)
expect_true(nrow(subset(get_transmissions(sir), source >= 0)) > 0)
expect_equal(nrow(el1), 1000L * 4L / 2L) # each undirected tie once

run(sir, ndays = 20, seed = 2)
expect_false(isTRUE(all.equal(el1, get_edgelist(contacts))))

run(sir, ndays = 20, seed = 1)
expect_equal(get_edgelist(contacts), el1)

# get_layer() gives the model's layer: edits through it show in `contacts`
expect_identical(get_edgelist(get_layer(sir, "contacts")), el1)
expect_error(get_layer(sir, "nope"))

# A fixed network: ties come from the edgelist; edits between runs persist ----
src <- c(0L, 1L, 2L)
tgt <- c(1L, 2L, 3L)
sis <- ModelSIS("cold", n = 4, prevalence = .25, transmission_prob = .5,
  recovery_rate = .5)
verbose_off(sis)
ring <- layer_network("ring", graph_edgelist(src, tgt, size = 4))
add_layer(sis, ring)

expect_true(has_edge(ring, 0, 1))
expect_true(has_edge(ring, 1, 0)) # undirected
expect_false(has_edge(ring, 0, 3))

expect_true(add_edge(ring, 0, 3))
run(sis, ndays = 5, seed = 1)
expect_true(has_edge(ring, 0, 3))
expect_equal(nrow(get_edgelist(ring)), 4L)

expect_true(rm_edge(ring, 0, 3))
expect_false(has_edge(ring, 0, 3))

# Network functions need a network layer, not a model
expect_error(get_edgelist(sis), "get_layer")

# Mixing layers ----------------------------------------------------------------
cm <- matrix(c(.9, .1, .1, .9), 2) * 10
seir <- ModelSEIR("flu", n = 1000, prevalence = .01,
  transmission_prob = .1, latent_period = 3, recovery_rate = .2)
add_layer(seir, layer_mixing("mixing", rep(0:1, each = 500), contact_matrix = cm))
mixing <- get_layer(seir, "mixing")
expect_inherits(mixing, "epiworld_layer_mixing")
expect_equal(get_contact_matrix(mixing), cm)

cm2 <- cm[2:1, ]
set_contact_matrix(mixing, cm2)
expect_equal(get_contact_matrix(get_layer(seir, "mixing")), cm2)

expect_error(get_contact_matrix(seir), "get_layer")
expect_error(get_edgelist(mixing), "network")

# Several layers in one draw: the layer column, weights and targets -----------
# A mixing layer with a contact matrix of zeros transmits nothing, so the model
# only transmits through the network (the mixing layer is layer 0, the network
# layer 1).
seir <- ModelSEIR(
  "flu",
  n                 = 1000,
  prevalence        = .05,
  transmission_prob = .5,
  latent_period     = 1,
  recovery_rate     = .1
)
add_layer(seir, layer_mixing("mixing", rep(0, 1000), contact_matrix = 0))
verbose_off(seir)
add_layer(seir, layer_network("household", graph_smallworld(k = 4, p = .1)))

run(seir, ndays = 20, seed = 1)
tr <- subset(get_transmissions(seir), source >= 0)
expect_true(nrow(tr) > 0)
expect_equal(unique(tr$layer), 1L)
expect_equal(unique(subset(get_transmissions(seir), source < 0)$layer), -1L)
expect_stdout(print(get_layer(seir, 1)), "household")

# Weight 0 turns the layer off
add_param(seir, "household weight", 0)
set_weight(get_layer(seir, "household"), "household weight")
run(seir, ndays = 20, seed = 1)
expect_equal(nrow(subset(get_transmissions(seir), source >= 0)), 0L)

# A layer that targets another virus does not spread this one
set_param(seir, "household weight", 1)
run(seir, ndays = 20, seed = 1)
expect_true(nrow(subset(get_transmissions(seir), source >= 0)) > 0)

other <- virus(
  "other",
  prevalence     = 0,
  as_proportion  = TRUE,
  prob_infecting = .5,
  recovery_rate  = .1
)
add_virus(seir, other)
add_virus_target(get_layer(seir, "household"), other)
run(seir, ndays = 20, seed = 1)
expect_equal(nrow(subset(get_transmissions(seir), source >= 0)), 0L)

# A contact rate of 0 silences the random mixing; the network still transmits
conn <- ModelSIR("flu", 1000, .05,
  transmission_prob = 1, recovery_rate = .1)
add_param(conn, "Contact rate", 0)
add_layer(
  conn,
  layer_mixing("everyone", groups = rep(0, 1000), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)
verbose_off(conn)
add_layer(conn, layer_network("contacts", graph_smallworld(k = 4, p = .1)))
run(conn, ndays = 10, seed = 1)
expect_equal(unique(subset(get_transmissions(conn), source >= 0)$layer), 1L)

# Mixing within groups: a contact probability per pair -------------------------
# 100 households of 10. With probability 1 the households are cliques: one
# carrier, transmitting for sure, infects its 9 housemates in a day, and nobody
# else.
sir <- ModelSIR("flu", n = 1000, prevalence = 0, transmission_prob = 1,
  recovery_rate = 0)
verbose_off(sir)
set_distribution_virus(get_virus(sir, 0), distribute_virus_to_set(0L))

hh <- layer_mixing(
  "households",
  groups       = (0:999) %/% 10,
  contact_prob = 1
)
expect_inherits(hh, "epiworld_layer_mixing")
add_layer(sir, hh)

infected <- function(m) which(get_agents_states(m) == "Infected") - 1L

run(sir, ndays = 1, seed = 1)
expect_equal(infected(sir), 0:9)
expect_equal(sort(members_of(hh, 1)), 10:19)
expect_equal(get_contact_prob(hh), rep(1, 100))

# Moves between runs are kept: the carrier joins household 1, then leaves
set_group(hh, 0, 1)
run(sir, ndays = 1, seed = 1)
expect_equal(infected(sir), c(0L, 10:19))
expect_equal(get_group(hh, 0), 1L)

set_group(hh, 0, NA)
run(sir, ndays = 1, seed = 1)
expect_equal(infected(sir), 0L)
expect_true(is.na(get_group(hh, 0)))

# A probability set from a global event (at the end of each day) lasts for
# that run only
set_group(hh, 0, 0)
add_globalevent(
  sir,
  globalevent_fun(
    function(m) set_contact_prob(get_layer(m, "households"), 0),
    "lockdown"
  )
)
run(sir, ndays = 1, seed = 1)
expect_equal(get_contact_prob(hh), rep(0, 100))
rm_globalevent(sir, "lockdown")
run(sir, ndays = 1, seed = 1)
expect_equal(get_contact_prob(hh), rep(1, 100))
expect_equal(infected(sir), 0:9)

# A traced layer logs its contacts: every transmission is in the log
expect_false(is_tracing(hh))
expect_error(get_tracing(hh), "not traced")
enable_tracing(hh, max_contacts = 20)
expect_true(is_tracing(hh))
run(sir, ndays = 1, seed = 1)
tr <- get_tracing(hh)
expect_equal(names(tr), c("source", "contact", "day"))
tx <- subset(get_transmissions(sir), source >= 0)
expect_equal(nrow(tx), 9L)
expect_true(all(paste(tx$source, tx$target) %in% paste(tr$source, tr$contact)))

# ...and keeps who shared a group: agent 0's nine housemates on day 1
gm <- groupmates_of(hh, 0)
expect_equal(names(gm), c("agent", "last_day", "n_days"))
expect_equal(gm$agent, 1:9)
expect_true(all(gm$last_day == 1L & gm$n_days == 1L))

# The history covers whole runs: a move between runs drops it
set_group(hh, 1, 1)
expect_error(groupmates_of(hh, 0), "no history")
set_group(hh, 1, 0)

# Moves during a run count from the next day. Two households of three: agent
# 2 moves to the second on day 3 and back on day 6; agent 4 leaves the layer
# on day 5 and returns on day 7.
trio <- ModelSIR("flu", n = 6, prevalence = 0, transmission_prob = .5,
  recovery_rate = .1)
verbose_off(trio)
homes <- layer_mixing("homes", groups = c(0, 0, 0, 1, 1, 1), contact_prob = 1)
enable_tracing(homes)
add_layer(trio, homes)
moves <- list(c(3, 2, 1), c(5, 4, NA), c(6, 2, 0), c(7, 4, 1))
for (mv in moves)
  add_globalevent(
    trio,
    globalevent_fun(
      local({
        mv <- mv
        function(m) set_group(get_layer(m, "homes"), mv[2], mv[3])
      }),
      name = paste("move", mv[2], "on day", mv[1]),
      day  = mv[1]
    )
  )
run(trio, ndays = 8, seed = 1)

# Agent 2: home on days 1-3 and 7-8, in the second household on days 4-6,
# where agent 4 was on days 4-5 (gone on 6-7, back on 8)
expect_equal(
  groupmates_of(homes, 2),
  data.frame(
    agent    = c(0L, 1L, 3L, 4L, 5L),
    last_day = c(8L, 8L, 6L, 5L, 6L),
    n_days   = c(5L, 5L, 3L, 2L, 3L)
  )
)
expect_equal(
  groupmates_of(homes, 4, from = 5, to = 8),
  data.frame(
    agent    = c(2L, 3L, 5L),
    last_day = c(5L, 8L, 8L),
    n_days   = c(1L, 2L, 2L)
  )
)
disable_tracing(hh)
expect_false(is_tracing(hh))

# Probability .5 per housemate is the same run as 4.5 contacts a day with the
# 9 housemates, from a diagonal contact matrix
same_run <- function(...) {
  m <- ModelSIR("flu", n = 1000, prevalence = .02, transmission_prob = .3,
    recovery_rate = .1)
  verbose_off(m)
  add_layer(m, layer_mixing("households", (0:999) %/% 10, ...))
  run(m, ndays = 30, seed = 2026)
  get_agents_states(m)
}
expect_equal(
  same_run(contact_prob = .5),
  same_run(contact_matrix = diag(4.5, 100))
)

expect_error(layer_mixing("x", 0:9), "exactly one")
expect_error(layer_mixing("x", c(0, .5), contact_prob = 1), "groups")
expect_error(layer_mixing("x", c(0, -1), contact_prob = 1), "groups")
expect_error(set_contact_prob(hh, 1.5))
expect_error(get_group(hh, 1000), "out of range")
expect_error(members_of(hh, 100), "out of range")

# A layer whose groups do not cover the model's agents fails when it runs
bad <- ModelSIR("flu", n = 1000, prevalence = .01, transmission_prob = 1,
  recovery_rate = 0)
verbose_off(bad)
add_layer(bad, layer_mixing("y", 0:9, contact_prob = 1))
expect_error(run(bad, ndays = 1), "agents")

# Groups drawn at random at every run: the sizes hold (agents left over are in
# no group), the same seed draws the same groups and another seed others, and
# a move between runs is an error (the next run would draw them again)
schools <- layer_mixing(
  "schools",
  groups       = grouping_random(c(300, 200)),
  contact_prob = .05
)
drawn <- ModelSIR("flu", n = 1000, prevalence = .01, transmission_prob = .2,
  recovery_rate = .1)
verbose_off(drawn)
add_layer(drawn, schools)
draw <- function(seed) {
  run(drawn, ndays = 10, seed = seed)
  sapply(0:999, get_group, layer = schools)
}
first <- draw(1)
expect_equal(as.vector(table(first, useNA = "always")), c(300L, 200L, 500L))
expect_identical(draw(1), first)
expect_false(identical(draw(2), first))
expect_error(set_group(schools, 0, 1), "drawn at every reset")

# Proportions are of the agents still without a group
halves <- layer_mixing(
  "halves",
  groups       = grouping_random(c(.5, .5), as_proportion = TRUE),
  contact_prob = .05
)
add_layer(drawn, halves)
run(drawn, ndays = 1, seed = 1)
expect_equal(lengths(list(members_of(halves, 0), members_of(halves, 1))),
  c(500L, 250L))
expect_error(grouping_random(-1), "non-negative")
expect_error(grouping_random(1.9), "whole numbers")

# Tools and viruses by group: a share of each group of a mixing layer, here a
# vaccine for 25% and 50% of two schools (100 + 300 agents) and the 10 initial
# cases in the first
schools <- ModelSIR("flu", n = 1000, prevalence = 0, transmission_prob = 0,
  recovery_rate = 0)
verbose_off(schools)
add_layer(
  schools,
  layer_mixing("schools", groups = rep(0:1, c(400, 600)), contact_prob = 0)
)
mmr <- tool(
  "MMR",
  prevalence               = 0,
  as_proportion            = TRUE,
  susceptibility_reduction = 1,
  transmission_reduction   = 0,
  recovery_enhancer        = 0,
  death_reduction          = 0
)
set_distribution_tool(mmr, distribute_tool_to_groups("schools", c(.25, .5)))
add_tool(schools, mmr)
set_distribution_virus(
  get_virus(schools, 0),
  distribute_virus_to_groups("schools", c(10, 0), as_proportion = FALSE)
)
run(schools, ndays = 1, seed = 1)
vaxed <- subset(get_hist_tool(schools), date == 0)
expect_equal(sum(vaxed$counts), 100L + 300L)
infected <- which(get_agents_states(schools) == "Infected")
expect_equal(length(infected), 10L)
expect_true(all(infected <= 400))
expect_error(distribute_tool_to_groups(1, .5))
expect_error(distribute_tool_to_groups("schools", 1.5, as_proportion = FALSE))

# A contact rate from a model parameter: a global event on day 10 that sets it
# to 0 stops new infections from day 11 on
crate <- ModelSIR("flu", n = 2000, prevalence = .01, transmission_prob = .1,
  recovery_rate = 1 / 14)
verbose_off(crate)
add_param(crate, "Contact rate", 5)
add_layer(
  crate,
  layer_mixing("everyone", groups = rep(0, 2000), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
)
add_globalevent(crate, globalevent_set_params("Contact rate", 0, day = 10))
run(crate, ndays = 30, seed = 2211)
expect_equal(max(subset(get_transmissions(crate), source >= 0)$date), 10)

# Infectious states: Exposed agents that never leave Exposed infect nobody,
# until every state is made infectious
seir <- ModelSEIR("flu", n = 2000, prevalence = .01, transmission_prob = .5,
  latent_period = 1e9, recovery_rate = 0)
verbose_off(seir)
expect_equal(get_infectious_states(seir), 2L)
add_layer(
  seir,
  layer_mixing("everyone", groups = rep(0, 2000), contact_matrix = 5)
)
run(seir, ndays = 5, seed = 1)
expect_equal(nrow(subset(get_transmissions(seir), source >= 0)), 0L)
set_infectious_states(seir, integer(0))
run(seir, ndays = 5, seed = 1)
expect_true(nrow(subset(get_transmissions(seir), source >= 0)) > 0L)

# Defunct functions say what replaces them -------------------------------------
expect_error(agents_smallworld(sir, 1000, 5, FALSE, .01), "layer_network")
expect_error(agents_from_edgelist(sir), "graph_edgelist")
expect_error(get_network(sir), "get_edgelist")
expect_error(ModelSIRCONN("x", 100, .1, 2, .5, .3), "set_contact_rate")
expect_error(ModelSEIRCONN(), "ModelSEIR\\(\\)")
expect_error(ModelSIRDCONN(), "ModelSIRD\\(\\)")
expect_error(ModelSEIRDCONN(), "ModelSEIRD\\(\\)")
expect_error(ModelSIRMixing(), "layer_mixing")
expect_error(ModelSEIRMixing(), "layer_mixing")
expect_error(ModelSEIRMixingQuarantine(), "ModelSEIRQuarantine")
expect_error(update_fun_susceptible(exclude = 1L), "set_infectious_states")
