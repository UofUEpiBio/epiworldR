# Schedules, persistence and parameter cycles (epiworld Phase 4c).

weekdays <- c(rep(TRUE, 5), FALSE, FALSE)
weekend <- function(day) (day - 1) %% 7 >= 5

# A work layer on weekdays only: no transmission through it on weekends,
# while the home layer goes on.
model <- ModelSIR(
  name              = "flu",
  n                 = 1000,
  prevalence        = 0.02,
  transmission_prob = 0.3,
  recovery_rate     = 0.1
)
verbose_off(model)
home <- layer_network("home", graph_smallworld(k = 4, p = 0.1))
work <- layer_network("work", graph_smallworld(k = 6, p = 0.1))
add_layer(model, home)
add_layer(model, work)

set <- withVisible(set_schedule(work, weekdays))
expect_identical(set$value, work)
expect_false(set$visible)
expect_true(has_schedule(work))
expect_false(has_schedule(home))
expect_true(is_scheduled(work, 5))
expect_false(is_scheduled(work, 6))
expect_true(is_scheduled(work, 8))
expect_error(is_scheduled(work, 0))

run(model, ndays = 40, seed = 11)
tx <- subset(get_transmissions(model), source >= 0)
expect_equal(sum(tx$layer == 1L & weekend(tx$date)), 0L)
expect_true(sum(tx$layer == 1L) > 10L)
expect_true(sum(tx$layer == 0L & weekend(tx$date)) > 0L)

clear_schedule(work)
expect_false(has_schedule(work))
expect_error(set_schedule(work, logical(0)), "logical vector")
expect_error(set_schedule(work, c(TRUE, NA)), "logical vector")

# Persistence: 50 trips of 60 riders on 10 routes of 5 trips. With
# persistence 0 every rider takes another trip of its route each day, so
# transit transmissions cross trips but never routes.
n <- 3000L
trip <- function(id) id %/% 60L
route <- function(t) t %/% 5L
model <- ModelSIR(
  name              = "flu",
  n                 = n,
  prevalence        = 0.01,
  transmission_prob = 0.3,
  recovery_rate     = 0.1
)
verbose_off(model)
add_layer(model, layer_network("home", graph_smallworld(k = 4, p = 0.05)))
transit <- layer_mixing("transit", groups = trip(0:(n - 1)), contact_prob = 0.05)
set_pools(transit, route(0:49))
set_persistence(transit, 0)
add_layer(model, transit)
expect_identical(get_pools(transit), route(0:49))
expect_identical(get_persistence(transit), 0)

run(model, ndays = 30, seed = 3)
tx <- subset(get_transmissions(model), source >= 0 & layer == 1L)
expect_true(nrow(tx) > 10L)
expect_true(all(route(trip(tx$source)) == route(trip(tx$target))))
expect_true(any(trip(tx$source) != trip(tx$target)))
expect_identical(get_usual_group(transit, 0), 0L)
expect_true(get_group(transit, 0) != 0L)

clear_persistence(transit)
expect_identical(get_persistence(transit), numeric(0))
expect_error(set_persistence(transit, 1.5))
expect_error(set_pools(transit, -1))

# A weekly contact rate, and the same run from a global event that sets the
# parameter at the end of each day (the value of the next day).
rates <- c(2, 2, 2, 2, 2, 5, 5)
public_model <- function() {
  m <- ModelSIR(
    name              = "flu",
    n                 = 1000,
    prevalence        = 0.01,
    transmission_prob = 0.05,
    recovery_rate     = 0.1
  )
  verbose_off(m)
  add_param(m, "Public rate", 2)
  public <- layer_mixing("public", groups = rep(0L, 1000), contact_matrix = matrix(1))
  set_contact_rate(public, "Public rate")
  add_layer(m, public)
  m
}

cycled <- public_model()
set <- withVisible(set_param_cycle(cycled, "Public rate", rates))
expect_identical(set$value, cycled)
expect_false(set$visible)
expect_true(has_param_cycle(cycled, "Public rate"))
expect_false(has_param_cycle(cycled, "Recovery rate"))

by_hand <- public_model()
add_globalevent(by_hand, globalevent_fun(function(model) {
  set_param(model, "Public rate", rates[today(model) %% 7 + 1])
}, "Cycle"))

run(cycled, ndays = 30, seed = 7)
run(by_hand, ndays = 30, seed = 7)
expect_true(nrow(get_transmissions(cycled)) > 50L)
expect_identical(get_transmissions(cycled), get_transmissions(by_hand))
expect_equal(get_param(cycled, "Public rate"), rates[(30 - 1) %% 7 + 1])

clear_param_cycle(cycled, "Public rate")
expect_false(has_param_cycle(cycled, "Public rate"))
expect_error(set_param_cycle(cycled, "Public rate", numeric(0)), "numeric vector")
expect_error(set_param_cycle(cycled, "No such parameter", 1), "does not exist")
