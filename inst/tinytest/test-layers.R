# Contact layers: network layers from generators and edgelists, the layer
# objects shared with the model, mixing layers (across and within groups), and
# the defunct functions.

# A model without layers has no contacts ---------------------------------------
sir <- ModelSIR("flu", n = 1000, prevalence = .05, transmission_rate = 1,
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
sis <- ModelSIS("cold", n = 4, prevalence = .25, transmission_rate = .5,
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
seir <- ModelSEIRMixing("flu", n = 1000, prevalence = .01,
  transmission_rate = .1, incubation_days = 3, recovery_rate = .2,
  contact_matrix = cm)
mixing <- get_layer(seir, "mixing")
expect_inherits(mixing, "epiworld_layer_mixing")
expect_equal(get_contact_matrix(mixing), cm)

cm2 <- cm[2:1, ]
set_contact_matrix(mixing, cm2)
expect_equal(get_contact_matrix(get_layer(seir, "mixing")), cm2)

expect_error(get_contact_matrix(seir), "get_layer")
expect_error(get_edgelist(mixing), "network")

# Several layers in one draw: the layer column, weights and targets -----------
# A mixing model with a contact matrix of zeros only transmits through the
# network added to it (the mixing layer is layer 0, the network layer 1).
seir <- ModelSEIRMixing(
  "flu",
  n                 = 1000,
  prevalence        = .05,
  transmission_rate = .5,
  incubation_days   = 1,
  recovery_rate     = .1,
  contact_matrix    = matrix(0, 1, 1)
)
add_entity(seir, entity("everyone", 1000, FALSE))
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

# The connected models ignore layers
conn <- ModelSIRCONN("flu", 1000, .05, contact_rate = 0,
  transmission_rate = 1, recovery_rate = .1)
verbose_off(conn)
add_layer(conn, layer_network("contacts", graph_smallworld(k = 4, p = .1)))
run(conn, ndays = 10, seed = 1)
expect_equal(nrow(subset(get_transmissions(conn), source >= 0)), 0L)

# Mixing within groups: a contact probability per pair -------------------------
# 100 households of 10. With probability 1 the households are cliques: one
# carrier, transmitting for sure, infects its 9 housemates in a day, and nobody
# else.
sir <- ModelSIR("flu", n = 1000, prevalence = 0, transmission_rate = 1,
  recovery_rate = 0)
verbose_off(sir)
set_distribution_virus(get_virus(sir, 0), distribute_virus_to_set(0L))

hh <- layer_mixing(
  "households",
  groups            = (0:999) %/% 10,
  infectious_states = 1,
  available_states  = 0:2,
  contact_prob      = 1
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

# Probability .5 per housemate is the same run as 4.5 contacts a day with the
# 9 housemates, from a diagonal contact matrix
same_run <- function(...) {
  m <- ModelSIR("flu", n = 1000, prevalence = .02, transmission_rate = .3,
    recovery_rate = .1)
  verbose_off(m)
  add_layer(m, layer_mixing("households", (0:999) %/% 10, 1, 0:2, ...))
  run(m, ndays = 30, seed = 2026)
  get_agents_states(m)
}
expect_equal(
  same_run(contact_prob = .5),
  same_run(contact_matrix = diag(4.5, 100))
)

expect_error(layer_mixing("x", 0:9, 1, 0:2), "exactly one")
expect_error(add_layer(sir, layer_mixing("y", 0:9, 1, 0:2, contact_prob = 1)),
  "agents")
expect_error(set_contact_prob(hh, 1.5))

# Defunct functions say what replaces them -------------------------------------
expect_error(agents_smallworld(sir, 1000, 5, FALSE, .01), "layer_network")
expect_error(agents_from_edgelist(sir), "graph_edgelist")
expect_error(get_network(sir), "get_edgelist")
