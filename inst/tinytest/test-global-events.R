# Test just this file: tinytest::run_test_file("inst/tinytest/test-global-events.R")

# Create a simple SIR model for testing
model <- ModelSIRCONN(
  name = "Test Model",
  n = 1000,
  prevalence = 0.01,
  contact_rate = 5,
  transmission_rate = 0.4,
  recovery_rate = 0.3
)

# Test adding a global event with globalevent_set_params
event1 <- globalevent_set_params(
  param = "Contact rate",
  value = 2,
  name = "Reduce Contact",
  day = 10
)

expect_inherits(event1, "epiworld_globalevent_set_param")
expect_inherits(event1, "epiworld_globalevent")

# Add the event to the model
expect_silent(add_globalevent(model, event1))

# Test adding a global event with globalevent_fun
test_function <- function(m) {
  # Do nothing, just a test
  invisible(NULL)
}

event2 <- globalevent_fun(
  fun = test_function,
  name = "Test Function Event",
  day = -99
)

expect_inherits(event2, "epiworld_globalevent_fun")
expect_inherits(event2, "epiworld_globalevent")

# Add the second event to the model
expect_silent(add_globalevent(model, event2))

# Test removing a global event by name
expect_silent(rm_globalevent(model, "Reduce Contact"))

# Test removing the second event
expect_silent(rm_globalevent(model, "Test Function Event"))

# Test with a tool-based global event
epitool <- tool(
  name = "Vaccine",
  prevalence = 0,
  as_proportion = FALSE,
  susceptibility_reduction = 0.9,
  transmission_reduction = 0.5,
  recovery_enhancer = 0.5,
  death_reduction = 0.9
)

event3 <- globalevent_tool(
  tool = epitool,
  prob = 0.2,
  name = "Vaccination Event",
  day = 20
)

expect_inherits(event3, "epiworld_globalevent_tool")
expect_inherits(event3, "epiworld_globalevent")

# Add the tool event
expect_silent(add_globalevent(model, event3))

# Remove it by name
expect_silent(rm_globalevent(model, "Vaccination Event"))

# Test that removing a non-existent event gives an error
# The C++ code throws an error if the event is not found
expect_error(rm_globalevent(model, "NonExistent Event"))

# Test deprecated 'action' parameter warning
expect_warning(
  add_globalevent(model, action = event1),
  "The argument `action` is deprecated"
)
