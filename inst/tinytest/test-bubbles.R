# Test just this file: tinytest::run_test_file("inst/tinytest/test-bubbles.R")

###############################################################################
# Households of 3 agents, in agent order.
###############################################################################
make_hh <- function(n, hh_size = 3L)
  rep(seq_len(ceiling(n / hh_size)), each = hh_size)[seq_len(n)]

###############################################################################
# Household flavor, group_size = 1: each household is its own bubble, so every
# transmission must stay within a household.
###############################################################################
model <- ModelSEIR(
  name = "Flu", prevalence = 0.05, transmission_rate = 0.3,
  incubation_days = 4.5, recovery_rate = 1 / 7
)
agents_smallworld(model, n = 600, k = 8, d = FALSE, p = 0.10)
hh <- make_hh(600)

expect_silent(
  bubbles(
    model, household_id = hh, flavor = "household",
    group_size = 1, transmission_factor = 1.0, start_day = 0
  )
)

expect_silent(run(model, ndays = 80, seed = 99))

tn  <- get_transmissions(model)
sec <- tn[tn$source >= 0, , drop = FALSE]
hh_of <- function(id) hh[id + 1L] # transmissions use 0-based agent ids

expect_true(nrow(sec) > 0)                                   # outbreak happened
expect_equal(sum(hh_of(sec$source) != hh_of(sec$target)), 0) # no cross-household

###############################################################################
# Control: without bubbles, transmission does cross households.
###############################################################################
model_ctl <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
agents_smallworld(model_ctl, n = 600, k = 8, d = FALSE, p = 0.10)
run(model_ctl, ndays = 80, seed = 99)
tn_ctl  <- get_transmissions(model_ctl)
sec_ctl <- tn_ctl[tn_ctl$source >= 0, , drop = FALSE]
expect_true(sum(hh_of(sec_ctl$source) != hh_of(sec_ctl$target)) > 0)

###############################################################################
# transmission_factor = 0 with a single all-encompassing bubble blocks
# everything; transmission_factor = 1 does not.
###############################################################################
n_hh <- length(unique(hh))
secondary <- function(tf) {
  m <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
  agents_smallworld(m, n = 600, k = 8, d = FALSE, p = 0.10)
  bubbles(m, household_id = hh, flavor = "household",
          group_size = n_hh, transmission_factor = tf, start_day = 0)
  run(m, ndays = 60, seed = 7)
  tn <- get_transmissions(m)
  sum(tn$source >= 0)
}
expect_equal(secondary(0.0), 0L)
expect_true(secondary(1.0) > 0)

###############################################################################
# Peer flavor and changing (rewiring) bubbles run without error.
###############################################################################
model_peer <- ModelSEIR("Flu", 0.05, 0.2, 4.5, 1 / 7)
agents_smallworld(model_peer, n = 600, k = 6, d = FALSE, p = 0.05)
expect_silent(
  bubbles(model_peer, household_id = hh, flavor = "peer", group_size = 1,
          transmission_factor = 0.5, start_day = 10, rewire_every = 7)
)
expect_silent(run(model_peer, ndays = 60, seed = 5))

###############################################################################
# Input validation.
###############################################################################
model_v <- ModelSEIR("Flu", 0.05, 0.2, 4.5, 1 / 7)
agents_smallworld(model_v, n = 600, k = 6, d = FALSE, p = 0.05)
expect_error(bubbles(model_v, household_id = hh[1:10], flavor = "household",
                     group_size = 2))                                   # wrong length
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = 2, transmission_factor = 2))          # bad factor
