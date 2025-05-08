# Test just this file: tinytest::run_test_file("inst/tinytest/test-virus.R")

# Check creating valid distribution function
expect_silent(distfun <- distribute_virus_randomly(
  prevalence = 2L,
  as_proportion = FALSE,
  agents_ids = integer(0)
))

expect_inherits(distfun, "epiworld_virus_distfun")
expect_length(class(distfun), 1)

# Check distributing virus to invalid agent IDs
expect_error(distribute_virus_randomly(
  prevalence = 2L,
  as_proportion = FALSE,
  agents_ids = c(-1, -2, -33)
), "ID must be a positive integer")