# Test just this file: tinytest::run_test_file("inst/tinytest/test-data.R")

# Test functions from R/data.R
# This test file validates all data access functions for epiworld models.

# Setup: Create a SEIRCONN model for testing -----------------------------------
model <- ModelSEIRCONN(
  name              = "Test Disease",
  n                 = 1000,
  prevalence        = 0.1,
  contact_rate      = 2.0,
  transmission_rate = 0.8,
  incubation_days   = 7.0,
  recovery_rate     = 0.3
)

verbose_off(model)
set.seed(1231)
run(model, ndays = 50)

# Test get_hist_total ----------------------------------------------------------
expect_silent(hist_total <- get_hist_total(model))
expect_inherits(hist_total, "epiworld_hist_total")
expect_inherits(hist_total, "epiworld_hist")
expect_inherits(hist_total, "data.frame")
expect_true(nrow(hist_total) > 0)
expect_true(all(c("date", "state", "counts") %in% names(hist_total)))
# Check that states attribute exists
expect_true(!is.null(attr(hist_total, "states")))

# Test get_today_total ---------------------------------------------------------
expect_silent(today_total <- get_today_total(model))
expect_true(is.numeric(today_total))
expect_true(length(today_total) > 0)
expect_true(!is.null(names(today_total)))

# Test get_hist_virus ----------------------------------------------------------
expect_silent(hist_virus <- get_hist_virus(model))
expect_inherits(hist_virus, "epiworld_hist_virus")
expect_inherits(hist_virus, "epiworld_hist")
expect_inherits(hist_virus, "data.frame")
expect_true(nrow(hist_virus) > 0)
expect_true(all(c("date", "state", "counts") %in% names(hist_virus)))
# Check that states attribute exists
expect_true(!is.null(attr(hist_virus, "states")))

# Test get_hist_tool -----------------------------------------------------------
expect_silent(hist_tool <- get_hist_tool(model))
expect_inherits(hist_tool, "epiworld_hist_tool")
expect_inherits(hist_tool, "epiworld_hist")
expect_inherits(hist_tool, "data.frame")
# hist_tool may be empty if no tools are added
expect_true(all(c("date", "state", "counts") %in% names(hist_tool)))

# Test get_transition_probability ----------------------------------------------
expect_silent(trans_prob <- get_transition_probability(model))
expect_true(is.matrix(trans_prob))
expect_true(nrow(trans_prob) == ncol(trans_prob))
expect_true(!is.null(dimnames(trans_prob)))
expect_true(!is.null(rownames(trans_prob)))
expect_true(!is.null(colnames(trans_prob)))

# Test get_reproductive_number -------------------------------------------------
expect_silent(repnum <- get_reproductive_number(model))
expect_inherits(repnum, "epiworld_repnum")
expect_inherits(repnum, "data.frame")
expect_true(nrow(repnum) > 0)
# Check expected columns
expect_true(all(c("source", "rt") %in% names(repnum)))

# Test plot.epiworld_repnum ----------------------------------------------------
expect_silent(repnum_plot <- plot(repnum, plot = FALSE))
expect_inherits(repnum_plot, "data.frame")
expect_true(all(c("virus_id", "virus", "date", "avg", "n", "sd", "lb", "ub") %in% names(repnum_plot)))

# Test that plot works with plot = TRUE (creates plot but no error)
expect_silent(plot(repnum, plot = TRUE))

# Test plot_reproductive_number ------------------------------------------------
expect_silent(plot_reproductive_number(model, plot = FALSE))

# Test get_hist_transition_matrix ----------------------------------------------
expect_silent(trans_matrix <- get_hist_transition_matrix(model))
expect_inherits(trans_matrix, "epiworld_hist_transition")
expect_inherits(trans_matrix, "data.frame")
expect_true(all(c("state_from", "state_to", "date", "counts") %in% names(trans_matrix)))
expect_true(!is.null(attr(trans_matrix, "states")))
expect_true(!is.null(attr(trans_matrix, "nsteps")))

# Test with skip_zeros = TRUE
expect_silent(trans_matrix_skip <- get_hist_transition_matrix(model, skip_zeros = TRUE))
expect_inherits(trans_matrix_skip, "epiworld_hist_transition")
expect_inherits(trans_matrix_skip, "data.frame")

# Test as.array.epiworld_hist_transition ---------------------------------------
expect_silent(trans_array <- as.array(trans_matrix))
expect_true(is.array(trans_array))
expect_equal(length(dim(trans_array)), 3)
expect_true(!is.null(dimnames(trans_array)))

# Test plot_incidence ----------------------------------------------------------
expect_silent(incidence_plot <- plot_incidence(model, plot = FALSE))
expect_inherits(incidence_plot, "data.frame")

# Test plot.epiworld_hist_transition -------------------------------------------
expect_silent(trans_plot <- plot(trans_matrix, plot = FALSE))
expect_inherits(trans_plot, "data.frame")

# Test get_transmissions -------------------------------------------------------
expect_silent(transmissions <- get_transmissions(model))
expect_inherits(transmissions, "epiworld_transmissions")
expect_inherits(transmissions, "data.frame")
expect_true(nrow(transmissions) > 0)
# Check expected columns according to documentation
expect_true(all(c("date", "source", "target", "virus_id", "virus", "source_exposure_date") %in% names(transmissions)))
# Check that seeded infections are included with source = -1
expect_true(any(transmissions$source == -1))

# Test get_generation_time -----------------------------------------------------
expect_silent(gentime <- get_generation_time(model))
expect_inherits(gentime, "epiworld_generation_time")
expect_inherits(gentime, "data.frame")
expect_true(nrow(gentime) > 0)
# Check expected columns according to documentation
expect_true(all(c("agent", "virus_id", "virus", "date", "gentime") %in% names(gentime)))
expect_true(!is.null(attr(gentime, "n_steps")))

# Test plot.epiworld_generation_time -------------------------------------------
expect_silent(gentime_plot <- plot(gentime, plot = FALSE))
expect_inherits(gentime_plot, "data.frame")

# Test plot_generation_time ----------------------------------------------------
expect_silent(plot_generation_time(model, plot = FALSE))

# Test get_active_cases --------------------------------------------------------
expect_silent(active_cases <- get_active_cases(model))
expect_inherits(active_cases, "data.frame")
expect_true(nrow(active_cases) > 0)
# Check expected columns according to documentation
expect_true(all(c("date", "virus_id", "virus", "active_cases") %in% names(active_cases)))

# Test get_outbreak_size -------------------------------------------------------
expect_silent(outbreak_size <- get_outbreak_size(model))
expect_inherits(outbreak_size, "data.frame")
expect_true(nrow(outbreak_size) > 0)
# Check expected columns according to documentation
expect_true(all(c("date", "virus_id", "virus", "outbreak_size") %in% names(outbreak_size)))

# Test with unrun model (should return empty results) --------------------------
model_unrun <- ModelSIRCONN(
  name = "Unrun Model",
  n = 100,
  prevalence = 0.01,
  contact_rate = 2,
  transmission_rate = 0.5,
  recovery_rate = 0.3
)

verbose_off(model_unrun)

# These should all return empty data but not error
expect_silent(hist_total_unrun <- get_hist_total(model_unrun))
expect_equal(nrow(hist_total_unrun), 0L)

expect_silent(hist_virus_unrun <- get_hist_virus(model_unrun))
expect_equal(nrow(hist_virus_unrun), 0L)

expect_silent(hist_tool_unrun <- get_hist_tool(model_unrun))
expect_equal(nrow(hist_tool_unrun), 0L)

expect_silent(repnum_unrun <- get_reproductive_number(model_unrun))
expect_equal(nrow(repnum_unrun), 0L)

expect_silent(trans_matrix_unrun <- get_hist_transition_matrix(model_unrun))
expect_equal(nrow(trans_matrix_unrun), 0L)

# Transition probability should still work even without running
expect_silent(trans_prob_unrun <- get_transition_probability(model_unrun))
expect_true(is.matrix(trans_prob_unrun))

# Test plot.epiworld_hist (generic plot for hist objects) ----------------------
expect_silent(plot(hist_total, plot = TRUE))
