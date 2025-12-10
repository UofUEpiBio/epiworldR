# Test just this file: tinytest::run_test_file("inst/tinytest/test-globalevents.R")

# =============================================================================
# Test 1: GlobalEvents - Event timing
# This test verifies that:
# 1. An event with a specific day is called only on that day
# 2. An event with day < 0 is called every day
#
# Note: Due to how the model initializes, global events start being called
# from day 1 (not day 0). Day 0 is the initial state before any updates.
# =============================================================================

# Create vectors to store the days when events are called
daily_event_days <- NULL
specific_event_days <- NULL

# Create a simple SIR model
model <- ModelSIRCONN(
  name              = "virus",
  n                 = 1000,
  prevalence        = 0.01,
  contact_rate      = 5,
  transmission_rate = 0.5,
  recovery_rate     = 0.3
)

verbose_off(model)

# Create a global event that runs every day (day = -99 means every day)
daily_event_fun <- function(m) {
  daily_event_days <<- c(daily_event_days, today(m))
}
daily_event <- globalevent_fun(daily_event_fun, "Daily event", day = -99)
expect_silent(add_globalevent(model, daily_event))

# Add global events that run on specific days (10, 20, 30)
# Note: the function is defined inside the loop because each event needs
# its own closure, but the behavior is identical
specific_event_fun <- function(m) {
  specific_event_days <<- c(specific_event_days, today(m))
}
for (day in c(10, 20, 30)) {
  specific_event <- globalevent_fun(
    specific_event_fun,
    paste0("Specific day event ", day),
    day = day
  )
  expect_silent(add_globalevent(model, specific_event))
}

# Run the model for 50 days
ndays <- 50
expect_silent(run(model, ndays = ndays, seed = 123))

# Verify daily event was called for all days (1 to ndays, inclusive)
expect_equal(
  length(daily_event_days),
  ndays,
  info = paste(
    "Daily event should be called", ndays, "times, but was called",
    length(daily_event_days), "times"
  )
)

# Verify daily event was called for each day (starting from 1, not 0)
expected_daily_days <- seq_len(ndays)
expect_equal(
  daily_event_days,
  expected_daily_days,
  info = "Daily event should be called on days 1 through ndays"
)

# Verify specific events were called exactly 3 times
expect_equal(
  length(specific_event_days),
  3,
  info = "Specific events should be called exactly 3 times"
)

# Verify specific events were called on correct days
expected_specific <- c(10, 20, 30)
expect_equal(
  specific_event_days,
  expected_specific,
  info = "Specific events should be called on days 10, 20, and 30"
)

# =============================================================================
# Test 2: GlobalEvents - Parameter modification
# This test verifies that:
# 1. A GlobalEvent can change the transmission rate to 0 on day 20
# 2. After the intervention, transmission should stop completely
# =============================================================================

# Test parameters
ndays <- 40
intervention_day <- 20
nexperiments <- 10

# Create a single SIR model
model2 <- ModelSIRCONN(
  name              = "virus",
  n                 = 1000,
  prevalence        = 0.05,
  contact_rate      = 8,
  transmission_rate = 0.7,
  recovery_rate     = 0.15
)

verbose_off(model2)

# Add a global event that sets transmission rate to 0 on the intervention day
stop_transmission <- globalevent_set_params(
  "Transmission rate",
  0,
  day = intervention_day
)
expect_silent(add_globalevent(model2, stop_transmission))

# Run multiple experiments
saver <- make_saver("transmission")
expect_silent(
  run_multiple(
    model2,
    ndays = ndays,
    nsims = nexperiments,
    seed = 1000,
    saver = saver,
    reset = TRUE,
    verbose = FALSE,
    nthreads = 1L
  )
)

# Retrieve results
results <- run_multiple_get_results(model2, nthreads = 1L)
trans_data <- results$transmission

# Count transmissions before and after intervention
# Transmissions on the intervention day are excluded since the global event
# runs after update_state() on that day
# Reset counts to 0 first to ensure we don't use stale values
total_transmissions_before <- 0
total_transmissions_after <- 0
if (nrow(trans_data) > 0) {
  total_transmissions_before <- sum(trans_data$date < intervention_day)
  total_transmissions_after <- sum(trans_data$date > intervention_day)
}

# We should have transmissions before the intervention
expect_true(
  total_transmissions_before > 0,
  info = "There should be transmissions before the intervention"
)

# After setting transmission to 0, there should be NO transmissions
expect_equal(
  total_transmissions_after,
  0,
  info = "There should be no transmissions after setting transmission rate to 0"
)

# =============================================================================
# Test 3: globalevent_tool - Adding a tool via global event
# This test verifies that tools can be distributed via global events
# =============================================================================

# Create a new model
model3 <- ModelSIRCONN(
  name              = "COVID-19",
  n                 = 10000,
  prevalence        = 0.01,
  contact_rate      = 5,
  transmission_rate = 0.4,
  recovery_rate     = 0.95
)

verbose_off(model3)

# Creating a tool (vaccine)
epitool <- tool(
  name = "Vaccine",
  prevalence = 0,
  as_proportion = FALSE,
  susceptibility_reduction = .9,
  transmission_reduction = .5,
  recovery_enhancer = .5,
  death_reduction = .9
)

# Adding a global event to distribute the tool on day 20
vaccine_event <- globalevent_tool(epitool, .2, day = 20)
expect_inherits(vaccine_event, "epiworld_globalevent_tool")
expect_inherits(vaccine_event, "epiworld_globalevent")
expect_silent(add_globalevent(model3, vaccine_event))

# Run the model
expect_silent(run(model3, ndays = 40, seed = 1912))

# Check that tools were distributed (via hist tool data)
tool_hist <- get_hist_tool(model3)
expect_true(
  nrow(tool_hist[tool_hist$date >= 20, ]) > 0,
  info = "Tool should have been distributed by the global event"
)

# =============================================================================
# Test 4: globalevent_fun - Custom function execution
# This test verifies that custom functions can be executed as global events
# =============================================================================

model4 <- ModelSIRCONN(
  name              = "SIR with Global Saver",
  n                 = 1000,
  prevalence        = 0.01,
  contact_rate      = 5,
  transmission_rate = 0.4,
  recovery_rate     = 0.3
)

verbose_off(model4)

# Create the object where we will store state counts
state_history <- NULL

# Function to record states
state_recorder <- function(m) {
  state_history <<- c(state_history, list(get_today_total(m)))
}

# Create and add the global event
recorder_event <- globalevent_fun(state_recorder, "State Recorder")
expect_inherits(recorder_event, "epiworld_globalevent_fun")
expect_inherits(recorder_event, "epiworld_globalevent")
expect_silent(add_globalevent(model4, recorder_event))

# Run the model
ndays4 <- 30
expect_silent(run(model4, ndays = ndays4, seed = 42))

# Verify we captured states for each day
expect_equal(
  length(state_history),
  ndays4,
  info = "State history should have one entry per simulation day"
)

# Each entry should be a named vector with state counts
expect_true(
  all(vapply(state_history, function(x) is.numeric(x) && !is.null(names(x)), logical(1))),
  info = "Each state history entry should be a named numeric vector"
)

# =============================================================================
# Test 5: rm_globalevent - Removing global events
# This test verifies that global events can be removed from a model
# =============================================================================

model5 <- ModelSIRCONN(
  name              = "Test removal",
  n                 = 1000,
  prevalence        = 0.01,
  contact_rate      = 5,
  transmission_rate = 0.4,
  recovery_rate     = 0.3
)

verbose_off(model5)

# Add two events
event1_calls <- 0
event2_calls <- 0

event1_fun <- function(m) {
  event1_calls <<- event1_calls + 1
}

event2_fun <- function(m) {
  event2_calls <<- event2_calls + 1
}

event1 <- globalevent_fun(event1_fun, "Event 1")
event2 <- globalevent_fun(event2_fun, "Event 2")

expect_silent(add_globalevent(model5, event1))
expect_silent(add_globalevent(model5, event2))

# Remove first event by name (rm_globalevent takes event name, not index)
expect_silent(rm_globalevent(model5, "Event 1"))

# Run model
expect_silent(run(model5, ndays = 10, seed = 123))

# Only event2 should have been called (since event1 was removed)
expect_equal(
  event1_calls,
  0,
  info = "Event 1 should not have been called after removal"
)

expect_equal(
  event2_calls,
  10,
  info = "Event 2 should have been called 10 times"
)

# =============================================================================
# Test 6: Print method for global events
# =============================================================================

# Test print for globalevent_fun
test_event <- globalevent_fun(function(m) NULL, "Test Print Event", day = 5)
expect_stdout(print(test_event))

# Test print for globalevent_set_param
param_event <- globalevent_set_params("Contact rate", 3.0, day = 15)
expect_stdout(print(param_event))

# Test print for globalevent_tool
tool_print <- tool(
  name = "PrintTest",
  prevalence = 0,
  as_proportion = FALSE,
  susceptibility_reduction = 0.5,
  transmission_reduction = 0.3,
  recovery_enhancer = 0.4,
  death_reduction = 0.5
)
tool_event <- globalevent_tool(tool_print, 0.5, day = 10)
expect_stdout(print(tool_event))

# =============================================================================
# Test 7: Deprecated function warnings
# =============================================================================

# Test that deprecated functions emit proper deprecation errors
expect_error(globalaction_tool(), "defunct")
expect_error(globalaction_tool_logit(), "defunct")
expect_error(globalaction_set_params(), "defunct")
expect_error(globalaction_fun(), "defunct")

# =============================================================================
# Test 8: add_globalevent with deprecated action parameter
# =============================================================================

model8 <- ModelSIRCONN(
  name              = "Test deprecated",
  n                 = 100,
  prevalence        = 0.1,
  contact_rate      = 2,
  transmission_rate = 0.3,
  recovery_rate     = 0.2
)

verbose_off(model8)

deprecated_event <- globalevent_fun(function(m) NULL, "Deprecated test")

# Using the deprecated 'action' parameter should produce a warning
expect_warning(
  add_globalevent(model8, action = deprecated_event),
  "deprecated"
)
