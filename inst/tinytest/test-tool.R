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

# Test tool probability setters ------------------------------------------------
tool_2 <- tool(
  name = "Test Tool",
  prevalence = 0.3,
  as_proportion = TRUE,
  susceptibility_reduction = 0.5,
  transmission_reduction = 0.0,
  recovery_enhancer = 0.0,
  death_reduction = 0.0
)

expect_silent(set_susceptibility_reduction(tool_2, 0.7))
expect_silent(set_transmission_reduction(tool_2, 0.4))
expect_silent(set_recovery_enhancer(tool_2, 0.3))
expect_silent(set_death_reduction(tool_2, 0.6))

# Test tool protection in epidemic simulation ----------------------------------
# Create model with highly infectious virus and protective tool for half agents
expect_silent(model_tool <- ModelSIRCONN(
  name = "Tool Protection Test",
  n = 100,
  prevalence = 0.1,
  contact_rate = 10,
  transmission_rate = 1.0,
  recovery_rate = 0.0
))

# Create a tool that fully protects agents (100% susceptibility reduction)
expect_silent(protective_tool <- tool(
  name = "Perfect Protection",
  prevalence = 0.5,
  as_proportion = TRUE,
  susceptibility_reduction = 1.0,
  transmission_reduction = 0.0,
  recovery_enhancer = 0.0,
  death_reduction = 0.0
))

expect_silent(add_tool(model_tool, protective_tool))
expect_silent(verbose_off(model_tool))
expect_silent(run(model_tool, ndays = 10, seed = 456))

hist_tool <- get_hist_total(model_tool)
final_data_tool <- hist_tool[hist_tool$date == max(hist_tool$date), ]
final_infected_tool <- final_data_tool[final_data_tool$state == "Infected", "counts"]
final_susceptible_tool <- final_data_tool[final_data_tool$state == "Susceptible", "counts"]

# With perfect protection for 50% of agents and highly infectious virus,
# approximately 50 agents should remain susceptible (protected)
# and approximately 50 should be infected (unprotected)
expect_true(final_susceptible_tool >= 40 && final_susceptible_tool <= 60,
            info = "Around half should remain protected")
expect_true(final_infected_tool >= 40 && final_infected_tool <= 60,
            info = "Around half should be infected")

# Test tool with pointer-based parameters --------------------------------------
expect_silent(model_ptr <- ModelSIRCONN(
  name = "Tool Ptr Test",
  n = 50,
  prevalence = 0.1,
  contact_rate = 5,
  transmission_rate = 0.8,
  recovery_rate = 0.2
))

# Add custom parameters to the model
expect_silent(add_param(model_ptr, "tool_efficacy", 0.85))
expect_silent(add_param(model_ptr, "recovery_boost", 0.4))
expect_silent(add_param(model_ptr, "death_protection", 0.7))

tool_ptr <- tool(
  name = "Pointer Tool",
  prevalence = 0.5,
  as_proportion = TRUE,
  susceptibility_reduction = 0.3,
  transmission_reduction = 0.0,
  recovery_enhancer = 0.0,
  death_reduction = 0.0
)

# Test pointer-based tool parameter setting
expect_silent(set_susceptibility_reduction_ptr(tool_ptr, model_ptr, "tool_efficacy"))
expect_silent(set_transmission_reduction_ptr(tool_ptr, model_ptr, "Transmission rate"))
expect_silent(set_recovery_enhancer_ptr(tool_ptr, model_ptr, "recovery_boost"))
expect_silent(set_death_reduction_ptr(tool_ptr, model_ptr, "death_protection"))

# Verify parameters can be retrieved
expect_equal(get_param(model_ptr, "tool_efficacy"), 0.85)
expect_equal(get_param(model_ptr, "recovery_boost"), 0.4)
expect_equal(get_param(model_ptr, "death_protection"), 0.7)

expect_silent(add_tool(model_ptr, tool_ptr))
expect_silent(verbose_off(model_ptr))
expect_silent(run(model_ptr, ndays = 5, seed = 890))
