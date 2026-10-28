# ModelSEIRQuarantine traces every layer it is given. A carrier (agent 0) in
# a ring network and in a household of four is detected the day it becomes
# infectious; nobody is infected (transmission probability 0), but its
# contacts that day are in the logs, so its two ring neighbors and its three
# housemates are quarantined, and nobody else.
n <- 40L

model <- ModelSEIRQuarantine(
  name                   = "flu",
  n                      = n,
  prevalence             = 0,
  transmission_prob      = 0,
  latent_period          = 1,
  recovery_rate          = 0,
  hospitalization_rate   = 0,
  hospitalization_period = 7,
  detection_delay        = 0,
  quarantine_period      = 14L,
  quarantine_willingness = 1,
  isolation_willingness  = 1,
  isolation_period       = 14L,
  contact_tracing_prob   = 1,
  contact_tracing_window = 4L
)
verbose_off(model)
set_distribution_virus(get_virus(model, 0), distribute_virus_to_set(0L))

ring <- layer_network(
  "ring",
  graph_edgelist(0:(n - 1L), c(1:(n - 1L), 0L), size = n)
)
households <- layer_mixing(
  "households",
  groups       = (0:(n - 1L)) %/% 4L,
  contact_prob = 1
)
add_layer(model, ring)
add_layer(model, households)

run(model, ndays = 3, seed = 1)

# The model turned tracing on in both layers
expect_true(is_tracing(ring))
expect_true(is_tracing(households))

states <- get_agents_states(model)
expect_equal(states[1], "Isolated")
expect_equal(
  which(states == "Quarantined Susceptible") - 1L,
  c(1L, 2L, 3L, n - 1L)
)

# The parameters have their new names
expect_equal(get_param(model, "Detection delay"), 0)
expect_equal(get_param(model, "Contact tracing window"), 4)

expect_error(ModelSEIRQuarantine(
  "flu", n, 0, 1.5, 1, 0, 0, 7, 0, 14L, 1, 1, 14L
), "transmission_prob")
