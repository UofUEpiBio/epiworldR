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
    group_size = 1, transmission_factor = 0.0, start_day = 0
  )
)

expect_silent(run(model, ndays = 80, seed = 99))

tn  <- get_transmissions(model)
sec <- tn[tn$source >= 0, , drop = FALSE]
hh_of <- function(id) hh[id + 1L] # transmissions use 0-based agent ids

expect_true(nrow(sec) > 0)                                   # outbreak happened
expect_equal(sum(hh_of(sec$source) != hh_of(sec$target)), 0) # no cross-household

###############################################################################
# Bubbles join households that are actually connected, so group_size = 2 must
# retain some cross-household transmission. (With tie-blind random pairing it
# would stay at ~0, making group_size indistinguishable from a lockdown.)
###############################################################################
model_g2 <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
agents_smallworld(model_g2, n = 600, k = 8, d = FALSE, p = 0.10)
bubbles(model_g2, household_id = hh, flavor = "household",
        group_size = 2, transmission_factor = 0.0, start_day = 0)
run(model_g2, ndays = 80, seed = 99)
tn_g2  <- get_transmissions(model_g2)
sec_g2 <- tn_g2[tn_g2$source >= 0, , drop = FALSE]
expect_true(sum(hh_of(sec_g2$source) != hh_of(sec_g2$target)) > 0)

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
# The transmission factor governs contact OUTSIDE the bubble, not inside it: 0
# cuts it entirely, 1 leaves the population untouched, and values in between are
# a soft contact reduction. Contact within the bubble is never altered.
#
# It also lives in the model, not in the intervention, so setting the parameter
# after bubbles() is what governs the run.
###############################################################################
transmissions_by_bubble <- function(tf, override = NULL) {
  m <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
  agents_smallworld(m, n = 600, k = 8, d = FALSE, p = 0.10)

  # group_size = 1: bubble == household, so cross-household transmission is
  # exactly the out-of-bubble transmission.
  bubbles(m, household_id = hh, flavor = "household",
          group_size = 1, transmission_factor = tf, start_day = 0)

  # bubbles() registers the factor as a model parameter...
  expect_equal(get_param(m, "Bubble transmission factor"), tf)

  # ...which the tool reads on every exposure.
  if (!is.null(override))
    set_param(m, "Bubble transmission factor", override)

  run(m, ndays = 60, seed = 7)
  tn  <- get_transmissions(m)
  sec <- tn[tn$source >= 0, , drop = FALSE]
  c(
    within = sum(hh_of(sec$source) == hh_of(sec$target)),
    outside = sum(hh_of(sec$source) != hh_of(sec$target))
  )
}

strict <- transmissions_by_bubble(0.0)
expect_equal(unname(strict[["outside"]]), 0L) # nothing crosses the bubble...
expect_true(strict[["within"]] > 0)           # ...but it still spreads inside

off <- transmissions_by_bubble(1.0)           # the policy imposes nothing
expect_true(off[["outside"]] > 0)

soft <- transmissions_by_bubble(0.25)         # a soft contact reduction
expect_true(soft[["outside"]] > 0)
expect_true(soft[["outside"]] < off[["outside"]])

# Changing the model parameter after bubbles() reproduces the tf = 1 run.
expect_equal(transmissions_by_bubble(0.0, override = 1.0), off)

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
# Peer bubbles must not swallow the population: max_households caps how far
# nominations chain. (Uncapped, the merges percolated into a single bubble that
# kept every contact, i.e. no restriction at all.)
###############################################################################
model_peer2 <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
agents_smallworld(model_peer2, n = 600, k = 8, d = FALSE, p = 0.10)
bubbles(model_peer2, household_id = hh, flavor = "peer", group_size = 2,
        transmission_factor = 0.0, start_day = 0, max_households = 2)
run(model_peer2, ndays = 60, seed = 21)

tn_p  <- get_transmissions(model_peer2)
sec_p <- tn_p[tn_p$source >= 0, , drop = FALSE]
# Some contact survives (households do pair up)...
expect_true(sum(hh_of(sec_p$source) != hh_of(sec_p$target)) > 0)
# ...but the policy still bites: far from every transmission crosses households,
# as it would if everyone shared one bubble.
expect_true(sum(hh_of(sec_p$source) == hh_of(sec_p$target)) > 0)

###############################################################################
# Replicates are independent: the partition belongs to the model, so the
# per-thread copies run_multiple() makes each keep their own. Running on two
# threads must give exactly the same replicates as running on one, and must not
# make the model single-threaded.
###############################################################################
model_mt <- ModelSEIR("Flu", 0.05, 0.3, 4.5, 1 / 7)
agents_smallworld(model_mt, n = 600, k = 8, d = FALSE, p = 0.10)
bubbles(model_mt, household_id = hh, flavor = "household", group_size = 2,
        transmission_factor = 0.0, start_day = 0)

expect_null(attr(model_mt, "single_threaded"))

outbreak_sizes <- function(m, nthreads) {
  saver <- make_saver("total_hist")
  # No warning about parallel execution: nothing is shared between copies.
  expect_silent(
    run_multiple(m, ndays = 40, nsims = 4, seed = 3311, saver = saver,
                 verbose = FALSE, nthreads = nthreads)
  )
  res <- run_multiple_get_results(m, nthreads = nthreads)$total_hist
  res <- res[res$date == 40 & res$state == "Recovered", ]
  res$counts[order(res$sim_num)]
}

serial_sizes   <- outbreak_sizes(model_mt, 1L)
parallel_sizes <- outbreak_sizes(model_mt, 2L)

expect_true(all(serial_sizes > 0))
expect_equal(parallel_sizes, serial_sizes)

###############################################################################
# Input validation.
###############################################################################
model_v <- ModelSEIR("Flu", 0.05, 0.2, 4.5, 1 / 7)
agents_smallworld(model_v, n = 600, k = 6, d = FALSE, p = 0.05)
expect_error(bubbles(model_v, household_id = hh[1:10], flavor = "household",
                     group_size = 2))                                   # wrong length
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = 2, transmission_factor = 2))          # bad factor
# Scalar arguments must be length-1, non-missing, whole numbers where integer.
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = 2, transmission_factor = c(0.5, 0.6))) # not scalar
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = NA))                                   # NA scalar
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = 2.5))                                  # non-whole int
expect_error(bubbles(model_v, household_id = hh, flavor = "peer",
                     group_size = 1, max_households = 1))                # cap too small
expect_error(bubbles(model_v, household_id = hh, flavor = "household",
                     group_size = 2, param_name = NA))                   # bad param name

###############################################################################
# A custom parameter name lets two policies be dialled independently.
###############################################################################
model_pn <- ModelSEIR("Flu", 0.05, 0.2, 4.5, 1 / 7)
agents_smallworld(model_pn, n = 600, k = 6, d = FALSE, p = 0.05)
bubbles(model_pn, household_id = hh, flavor = "household", group_size = 2,
        transmission_factor = 0.3, name = "Bubble A",
        param_name = "Bubble A factor")
expect_equal(get_param(model_pn, "Bubble A factor"), 0.3)
expect_silent(set_param(model_pn, "Bubble A factor", 0.6))
expect_equal(get_param(model_pn, "Bubble A factor"), 0.6)
expect_silent(run(model_pn, ndays = 30, seed = 12))
