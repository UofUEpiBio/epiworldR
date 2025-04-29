# Test just this file: tinytest::run_test_file("inst/tinytest/test-diffnet.R")

# Create DiffNet Model ---------------------------------------------------------
good_name <- "A Virus"
good_prevalence <- .01
good_prob_adopt <- .1
good_normalize_exposure <- TRUE
good_data <- matrix(0, nrow = 2, ncol = 2)
good_data_cols <- 1L:ncol(good_data)
good_params <- vector("double")

expect_silent(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  normalize_exposure = good_normalize_exposure,
  data = good_data,
  data_cols = good_data_cols,
  params = good_params
))

# Check model initialization
expect_inherits(diffnet_0, "epiworld_diffnet")
expect_inherits(diffnet_0, "epiworld_model")
expect_length(class(diffnet_0), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name <- NA
bad_prevalence <- NA
bad_prob_adopt <- NA
bad_normalize_exposure <- NA
bad_data <- NA
bad_data_cols <- NA
bad_params <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_bool <- "must be a boolean"
expected_error_msg_na <- "must not be NA"
expected_error_msg_any_na <- "must not contain NA values"
expected_error_msg_numvector <- "must be a numeric vector"

expect_error(diffnet_0 <- ModelDiffNet(
  name = bad_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt
), expected_error_msg_str)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = bad_prevalence,
  prob_adopt = good_prob_adopt
), expected_error_msg_na)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = bad_prob_adopt
), expected_error_msg_na)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  normalize_exposure = "string_bool"
), expected_error_msg_bool)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  normalize_exposure = bad_normalize_exposure
), expected_error_msg_na)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  data = bad_data
), expected_error_msg_na)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  data_cols = bad_data_cols
), expected_error_msg_numvector)

expect_error(diffnet_0 <- ModelDiffNet(
  name = good_name,
  prevalence = good_prevalence,
  prob_adopt = good_prob_adopt,
  params = bad_params
), expected_error_msg_numvector)


