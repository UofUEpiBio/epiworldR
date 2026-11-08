# Test just this file: tinytest::run_test_file("inst/tinytest/test-virus.R")

# A virus carries its natural history -----------------------------------------
sir_virus <- function(name = "Flu") {
  v <- virus(name)
  add_stage(v, "Susceptible")
  add_stage(v, "Infected", infectious = TRUE)
  add_stage(v, "Recovered")
  add_infection(v, "Susceptible", "Infected", rate = "Transmission probability")
  add_rate(v, "Infected", "Recovered", rate = "Recovery rate", effect = "recovery")
  v
}

flu <- sir_virus()
expect_inherits(flu, "epiworld_virus")
expect_equal(get_name_virus(flu), "Flu")
expect_silent(set_name_virus(flu, "Influenza"))
expect_equal(get_name_virus(flu), "Influenza")

expect_identical(
  get_stages(flu),
  data.frame(
    stage      = c("Susceptible", "Infected", "Recovered"),
    carrier    = c(FALSE, TRUE, FALSE),
    infectious = c(FALSE, TRUE, FALSE),
    contact    = c(1, 1, 1)
  )
)

# Stages are named by their labels, and checked
expect_error(add_rate(flu, "Infected", "Dead", rate = .1), "no stage \"Dead\"")
expect_error(add_rate(flu, "Infected", "Recovered", rate = NA), "-rate- must be")
expect_error(add_rate(flu, "Infected", "Recovered", rate = .1, effect = "x"), "effect")
expect_error(add_stage(flu, "X", infectious = TRUE, carrier = FALSE), "carrier")
expect_error(add_after(flu, "Infected", "Recovered", 3), "duration_")
expect_error(virus("Flu", 0.01, TRUE, 0.5), "takes only a name")
expect_error(set_prob_infecting(flu, .1), "defunct")
expect_error(virus_set_state(flu, 1, 2, 3), "defunct")

# Built in R, it is the same model as ModelSIR() ------------------------------
# The same parameters, stages and transitions give the same history, seed by
# seed.
built_in <- ModelSIR(
  name              = "Flu",
  n                 = 2000,
  prevalence        = 0.01,
  transmission_prob = 0.2,
  recovery_rate     = 0.1
)

from_scratch <- Model(2000)
add_param(from_scratch, "Recovery rate", 0.1)
add_param(from_scratch, "Transmission probability", 0.2)
add_virus(from_scratch, sir_virus(), assign_share(0.01))

for (m in list(built_in, from_scratch)) {
  add_layer(m, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
  verbose_off(m)
  run(m, ndays = 60, seed = 2231)
}

expect_identical(get_hist_total(built_in), get_hist_total(from_scratch))
expect_identical(get_stages(built_in), get_stages(from_scratch))
expect_true(
  subset(get_hist_total(from_scratch), date == 60 & stage == "Recovered")$counts > 100
)

# The SEIR's latent period is a duration
seir_virus <- virus("Flu")
add_stage(seir_virus, "Susceptible")
add_stage(seir_virus, "Exposed", carrier = TRUE)
add_stage(seir_virus, "Infected", infectious = TRUE)
add_stage(seir_virus, "Removed")
add_infection(seir_virus, "Susceptible", "Exposed", rate = "Transmission probability")
add_after(seir_virus, "Exposed", "Infected", duration_geometric("Latent period"))
add_rate(
  seir_virus, "Infected", "Removed", rate = "Recovery rate", effect = "recovery"
)

built_in <- ModelSEIR("Flu", 2000, 0.01, 0.2, 4, 0.1)
from_scratch <- Model(2000)
add_param(from_scratch, "Transmission probability", 0.2)
add_param(from_scratch, "Latent period", 4)
add_param(from_scratch, "Recovery rate", 0.1)
add_virus(from_scratch, seir_virus, assign_share(0.01))

for (m in list(built_in, from_scratch)) {
  add_layer(m, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
  verbose_off(m)
  run(m, ndays = 60, seed = 77)
}
expect_identical(get_hist_total(built_in), get_hist_total(from_scratch))

# A fixed duration: the seeds (in Exposed on day 0, with no transmission) all
# become infectious on day 4
fixed <- virus("Fixed")
add_stage(fixed, "Susceptible")
add_stage(fixed, "Exposed", carrier = TRUE)
add_stage(fixed, "Infected", infectious = TRUE)
add_infection(fixed, "Susceptible", "Exposed", rate = 0)
add_after(fixed, "Exposed", "Infected", duration_fixed(4))

m <- Model(100)
add_virus(m, fixed, assign_count(10))
verbose_off(m)
run(m, ndays = 6, seed = 1)
h <- subset(get_hist_total(m), stage == "Infected")
expect_equal(h$counts, c(0, 0, 0, 0, 10, 10, 10))

# Transitions by name: the outputs count each, and fire() applies one -------
vax <- sir_virus("Measles")
add_external(vax, "Susceptible", "Recovered", "Vaccination")

m <- Model(100)
add_param(m, "Recovery rate", 0.1)
add_param(m, "Transmission probability", 0)
add_virus(m, vax, assign_agents(0:4))
verbose_off(m)
run(m, ndays = 5, seed = 1)

expect_identical(
  unique(get_hist_transitions(m)$transition),
  c("Susceptible -> Infected", "Infected -> Recovered", "Vaccination")
)

expect_true(fire(m, "Vaccination", 50))
expect_false(fire(m, "Vaccination", 0)) # Agent 0 was infected: not susceptible
expect_error(fire(m, "Vaccinate", 50), "Vaccinate")

# Strains of the model's virus -------------------------------------------------
# A copy with its own name reads other parameters (use_param()): this one
# does not transmit, so only its seeds ever carry it.
model <- ModelSIR("A virus", 1000, 0.01, 0.3, 0.1)
add_layer(model, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
verbose_off(model)

strain <- clone_epi(get_virus(model, 0))
set_name_virus(strain, "B")
add_param(model, "B transmission", 0)
expect_silent(use_param(strain, "Transmission probability", "B transmission"))
add_virus(model, strain, assign_count(20))
run(model, ndays = 30, seed = 12)

hv <- get_hist_virus(model)
expect_identical(unique(hv$virus), c("A virus", "B"))
carried_b <- subset(hv, virus == "B" & stage != "Susceptible")
expect_true(all(tapply(carried_b$counts, carried_b$date, sum) <= 20))
# The virus history counts the agents that carry each virus (not the
# recovered): the model's virus spreads
carried_a <- subset(hv, virus == "A virus" & stage == "Infected")
expect_true(max(carried_a$counts) > 20)

# The natural history of a virus in a model cannot change
expect_error(add_stage(get_virus(model, 0), "Other"), "is in a model")

# Values ----------------------------------------------------------------------
# 1 / parameter, 1 - parameter, an R function, and a logit of the agents' data
expect_inherits(value_inverse("Infectious period"), "epiworld_value")
expect_inherits(value_complement("Efficacy"), "epiworld_value")

v <- virus("Values")
add_stage(v, "Susceptible")
add_stage(v, "Infected", infectious = TRUE)
add_stage(v, "Recovered")
add_infection(v, "Susceptible", "Infected", rate = value_complement("Efficacy"))
add_rate(v, "Infected", "Recovered", rate = value_inverse("Infectious period"))

m <- Model(500)
add_param(m, "Efficacy", 1) # 1 - 1: nobody is infected
add_param(m, "Infectious period", 1) # 1 / 1: everyone recovers the next day
add_layer(m, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
add_virus(m, v, assign_count(50))
verbose_off(m)
run(m, ndays = 3, seed = 1)
expect_equal(subset(get_hist_total(m), date == 3)$counts, c(450, 0, 50))

calls <- 0
v <- virus("R function")
add_stage(v, "Susceptible")
add_stage(v, "Infected", infectious = TRUE)
add_stage(v, "Recovered")
add_infection(v, "Susceptible", "Infected", rate = 0)
add_rate(
  v, "Infected", "Recovered",
  rate = value_fun(function(agent, model) {
    calls <<- calls + 1
    1
  })
)
expect_true(isTRUE(attr(v, "single_threaded")))
m <- Model(100)
add_virus(m, v, assign_count(10))
expect_true(isTRUE(attr(m, "single_threaded")))
verbose_off(m)
run(m, ndays = 2, seed = 1)
expect_equal(calls, 10) # Once per infected agent, on the day it recovers
expect_equal(subset(get_hist_total(m), date == 2)$counts, c(90, 0, 10))

m <- Model(100)
set_agents_data(m, matrix(rep(c(-100, 100), 50), ncol = 1))
v <- virus("Logit")
add_stage(v, "Susceptible")
add_stage(v, "Infected", infectious = TRUE)
add_stage(v, "Recovered")
add_infection(v, "Susceptible", "Infected", rate = 0)
add_rate(v, "Infected", "Recovered", rate = rate_fun_logit(0L, 1, m))
add_virus(m, v, assign_agents(0:9))
verbose_off(m)
run(m, ndays = 1, seed = 1)
# Agents with x = 100 recover (p ~ 1), those with x = -100 do not (p ~ 0)
expect_equal(get_agents_stages(m)[1:10], rep(c("Infected", "Recovered"), 5))
expect_error(rate_fun_logit(3L, 1, m), "out of range")

# Contact weights --------------------------------------------------------------
# A stage with contact weight 0 is out of every layer: its agents do not
# transmit
v <- sir_virus("Away")
set_contact(v, "Infected", 0)
expect_equal(get_stages(v)$contact, c(1, 0, 1))

m <- Model(1000)
add_param(m, "Recovery rate", 0)
add_param(m, "Transmission probability", 1)
add_layer(m, layer_network("contacts", graph_smallworld(k = 6, p = .05)))
add_virus(m, v, assign_count(10))
verbose_off(m)
run(m, ndays = 10, seed = 1)
expect_equal(subset(get_hist_total(m), date == 10 & stage == "Infected")$counts, 10)
