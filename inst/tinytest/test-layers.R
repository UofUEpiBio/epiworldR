# Contact layers: network layers from generators and edgelists, the layer
# objects shared with the model, mixing layers, and the defunct functions.

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

# Defunct functions say what replaces them -------------------------------------
expect_error(agents_smallworld(sir, 1000, 5, FALSE, .01), "layer_network")
expect_error(agents_from_edgelist(sir), "graph_edgelist")
expect_error(get_network(sir), "get_edgelist")
