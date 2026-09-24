# Test just this file: tinytest::run_test_file("inst/tinytest/test-hist-tool-names.R")

# get_hist_tool() should label rows with the tool names, not the virus names.

# One virus, one tool ----------------------------------------------------------
abm <- ModelSIRCONN(
  "myvirus", n = 10, prevalence = 0, contact_rate = 0,
  transmission_rate = 0, recovery_rate = 0
)

t1 <- tool("mytool", 0, FALSE, 1, 0, 0, 0)
set_distribution_tool(t1, distribute_tool_to_set(c(2L, 5L)))
add_tool(abm, t1)
run(abm, 0, seed = 1)

hist <- get_hist_tool(abm)
expect_true(nrow(hist) > 0L)
expect_equal(unique(hist$tool), "mytool")
expect_equal(unique(hist$tool_id), 0L)

# One virus, two tools ---------------------------------------------------------
abm2 <- ModelSIRCONN(
  "myvirus", n = 10, prevalence = 0, contact_rate = 0,
  transmission_rate = 0, recovery_rate = 0
)

t_a <- tool("tool_a", 0, FALSE, 1, 0, 0, 0)
set_distribution_tool(t_a, distribute_tool_to_set(c(2L, 5L)))
add_tool(abm2, t_a)

t_b <- tool("tool_b", 0, FALSE, 1, 0, 0, 0)
set_distribution_tool(t_b, distribute_tool_to_set(c(3L, 7L, 8L)))
add_tool(abm2, t_b)

run(abm2, 0, seed = 1)

hist2 <- get_hist_tool(abm2)
expect_equal(sort(unique(hist2$tool_id)), c(0L, 1L))
expect_equal(
  unique(hist2[hist2$tool_id == 0L, "tool"]), "tool_a"
)
expect_equal(
  unique(hist2[hist2$tool_id == 1L, "tool"]), "tool_b"
)
expect_false(any(hist2$tool == "myvirus"))
