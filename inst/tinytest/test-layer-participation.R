# Two contact layers, in order: home (0), work (1).
participation_model <- function() {
  model <- ModelSIR(
    name              = "flu",
    n                 = 1000,
    prevalence        = 0.02,
    transmission_prob = 0.5,
    recovery_rate     = 0.2
  )
  verbose_off(model)
  add_layer(model, layer_network("home", graph_smallworld(k = 4, p = 0.1)))
  add_layer(model, layer_network("work", graph_smallworld(k = 6, p = 0.2)))
  model
}

# Tool targets are names, deduplicated, and empty means every layer.
mask <- tool("Mask", 1, 0, 0, 0)
expect_identical(get_layer_targets_tool(mask), character(0))
added <- withVisible(add_layer_target_tool(mask, "work"))
expect_identical(added$value, mask)
expect_false(added$visible)
add_layer_target_tool(mask, "work")
expect_identical(get_layer_targets_tool(mask), "work")
set <- withVisible(set_layer_targets_tool(mask, c("home", "work", "home")))
expect_identical(set$value, mask)
expect_false(set$visible)
expect_identical(get_layer_targets_tool(mask), c("home", "work"))
set_layer_targets_tool(mask, NULL)
expect_identical(get_layer_targets_tool(mask), character(0))
set_layer_targets_tool(mask, character(0))
expect_identical(get_layer_targets_tool(mask), character(0))
add_layer_target_tool(mask, "work")
cleared <- withVisible(clear_layer_targets_tool(mask))
expect_identical(cleared$value, mask)
expect_false(cleared$visible)
expect_identical(get_layer_targets_tool(mask), character(0))

expect_error(add_layer_target_tool(mask, 1), "layer name")
expect_error(add_layer_target_tool(mask, c("home", "work")), "layer name")
expect_error(add_layer_target_tool(mask, NA_character_), "layer name")
expect_error(set_layer_targets_tool(mask, 1:2), "layer names")
expect_error(set_layer_targets_tool(mask, NA_character_), "layer names")
expect_error(get_layer_targets_tool(NULL), "epiworld_tool")
expect_error(clear_layer_targets_tool(NULL), "epiworld_tool")

# A mask at work blocks work only, in both transmission modes. Without layer
# targets it blocks all contacts. The initial infections have source < 0.
for (mode in c("pull", "push")) {
  model <- participation_model()
  set_transmission_mode(model, mode)
  mask <- tool("Mask", 1, 0, 0, 0)
  add_layer_target_tool(mask, "work")
  add_tool(model, mask, assign_share(1))
  expect_identical(get_layer_targets_tool(get_tool(model, 0)), "work")
  run(model, ndays = 30, seed = 99)
  tx <- subset(get_transmissions(model), source >= 0)
  expect_equal(sum(tx$layer == 1L), 0L)
  expect_true(sum(tx$layer == 0L) > 50L)

  # Change the model's copy; the original tool was copied by add_tool().
  clear_layer_targets_tool(get_tool(model, 0))
  run(model, ndays = 30, seed = 99)
  expect_equal(sum(get_transmissions(model)$source >= 0), 0L)
}

# SIR has one status, "Free" (everyone's): excluding it from work takes
# everyone out of work, while home still transmits.
model <- participation_model()
work <- get_layer(model, "work")
expect_identical(get_excluded_statuses(work), character(0))
excluded <- withVisible(exclude_status(work, "Free"))
expect_identical(excluded$value, work)
expect_false(excluded$visible)
expect_identical(get_excluded_statuses(work), "Free")
run(model, ndays = 30, seed = 99)
tx <- subset(get_transmissions(model), source >= 0)
expect_equal(sum(tx$layer == 1L), 0L)
expect_true(sum(tx$layer == 0L) > 0L)
expect_identical(get_excluded_statuses(work), "Free")

# Closing between runs persists across runs and logs no contacts. Reopening
# restores transmission through work.
model <- participation_model()
work <- get_layer(model, "work")
enable_tracing(work)
expect_false(is_layer_closed(work))
closed <- withVisible(close_layer(work))
expect_identical(closed$value, work)
expect_false(closed$visible)
expect_true(is_layer_closed(work))
for (seed in c(99, 100)) {
  run(model, ndays = 30, seed = seed)
  tx <- subset(get_transmissions(model), source >= 0)
  expect_equal(sum(tx$layer == 1L), 0L)
  expect_true(sum(tx$layer == 0L) > 0L)
  expect_equal(nrow(get_tracing(work)), 0L)
  expect_true(is_layer_closed(work))
}
opened <- withVisible(reopen_layer(work))
expect_identical(opened$value, work)
expect_false(opened$visible)
expect_false(is_layer_closed(work))
run(model, ndays = 30, seed = 99)
tx <- subset(get_transmissions(model), source >= 0)
expect_true(sum(tx$layer == 1L) > 0L)

# Names are resolved when the run starts, not when targets are set.
model <- participation_model()
mask <- tool("Mask", 1, 0, 0, 0)
add_layer_target_tool(mask, "missing")
add_tool(model, mask, assign_share(1))
expect_error(run(model, ndays = 30, seed = 99), 'targets the layer "missing"')

model <- participation_model()
exclude_status(get_layer(model, "work"), "missing")
expect_error(run(model, ndays = 30, seed = 99), 'excludes the status "missing"')

# Participation functions accept mixing layers too, and require layer objects.
mixing <- layer_mixing("mixing", groups = rep(0, 1000), contact_prob = 1)
exclude_status(mixing, "Free")
expect_identical(get_excluded_statuses(mixing), "Free")
close_layer(mixing)
expect_true(is_layer_closed(mixing))
reopen_layer(mixing)
expect_false(is_layer_closed(mixing))
expect_error(exclude_status(model, "Free"), "get_layer")
expect_error(exclude_status(mixing, NA_character_), "status name")
expect_error(get_excluded_statuses(model), "get_layer")
expect_error(close_layer(model), "get_layer")
expect_error(reopen_layer(model), "get_layer")
expect_error(is_layer_closed(model), "get_layer")
