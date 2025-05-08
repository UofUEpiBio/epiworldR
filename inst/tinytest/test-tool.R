# Test just this file: tinytest::run_test_file("inst/tinytest/test-tool.R")

# Create tool ------------------------------------------------------------------
tool_name <- "Vaccine"

expect_silent(tool_0 <- tool(
  name = tool_name,
  prevalence = 0.5,
  as_proportion = TRUE,
  susceptibility_reduction = .9,
  transmission_reduction = .5,
  recovery_enhancer = .5,
  death_reduction = .9
))

# Check tool initialization
expect_inherits(tool_0, "epiworld_tool")
expect_length(class(tool_0), 1)

# Check tool helper functions --------------------------------------------------
expect_stdout(print(tool_0))

expect_equal(get_name_tool(tool_0), tool_name)

new_tool_name <- "New Vaccine"
expect_silent(tool_renamed <- set_name_tool(tool_0, new_tool_name))
expect_equal(get_name_tool(tool_0), new_tool_name)
expect_equal(get_name_tool(tool_renamed), new_tool_name)

# Check adding tool to model ---------------------------------------------------
sir_0 <- ModelSIR(name = "COVID-19", prevalence = .01,
                  transmission_rate = .9, recovery_rate = .3)
expect_silent(sir_with_tool <- add_tool(sir_0, tool_0))
expect_inherits(sir_with_tool, "epiworld_model")

# Check adding tool with proportion parameter (deprecated)
expect_warning(add_tool(sir_0, tool_0, 3), "deprecated")

# Check removing tools
expect_error(rm_tool(sir_0, 2), "only 2 tools")
expect_silent(rm_tool(sir_0, 1))
expect_silent(rm_tool(sir_0, 0))
expect_error(rm_tool(sir_0, 0), "only 0 tools")

# Check distributing tool to invalid agent IDs
expect_error(distribute_tool_randomly(
  prevalence = 2L,
  as_proportion = FALSE,
  agents_ids = c(-1, -2, -33)
), "ID must be a positive integer")

# Check initialization with missing prevalence parameter (deprecated) ----------
expect_warning(tool_1 <- tool(
  name = "Vaccine",
  as_proportion = TRUE,
  susceptibility_reduction = .9,
  transmission_reduction = .5,
  recovery_enhancer = .5,
  death_reduction = .9
))

expect_true(attributes(tool_1)$uses_deprecated)