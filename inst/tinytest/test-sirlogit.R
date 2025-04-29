# Test just this file: tinytest::run_test_file("inst/tinytest/test-sirlogit.R")

# Create SIRLogit Model --------------------------------------------------------
good_name <- "A Virus"
good_data <- cbind(
  Intercept = 1,
  Female    = sample.int(2, 100000, replace = TRUE) - 1
)
good_coef_infect        <- c(.1, -2, 2)
good_coef_recover       <- rnorm(2)
good_coef_infect_cols   <- 1L:ncol(good_data)
good_coef_recover_cols  <- 1L:ncol(good_data)
good_prob_infection     <- .8
good_recovery_rate      <- .3
good_prevalence         <- .01

sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
)

# Check model initialization
expect_inherits(sirlogit_0, "epiworld_sir")
expect_inherits(sirlogit_0, "epiworld_model")
expect_length(class(sirlogit_0), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name                <- 10
bad_numeric_input       <- "not a number"
bad_data                <- c(1, 0, NA)
bad_coef_infect         <- NA
bad_coef_recover        <- NA
bad_coef_infect_cols    <- NA
bad_coef_recover_cols   <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_any_na <- "must not contain NA values"
expected_error_msg_numvector <- "must be a numeric vector"
expected_error_msg_double    <- "must be a double"

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = bad_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_str)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = bad_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_any_na)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = bad_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_numvector)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = bad_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_numvector)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = bad_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_numvector)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = bad_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_numvector)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = bad_numeric_input,
  recovery_rate     = good_recovery_rate,
  prevalence        = good_prevalence
), expected_error_msg_double)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = bad_numeric_input,
  prevalence        = good_prevalence
), expected_error_msg_double)

expect_error(sirlogit_0 <- ModelSIRLogit(
  vname             = good_name,
  data              = good_data,
  coefs_infect      = good_coef_infect,
  coefs_recover     = good_coef_recover,
  coef_infect_cols  = good_coef_infect_cols,
  coef_recover_cols = good_coef_recover_cols,
  prob_infection    = good_prob_infection,
  recovery_rate     = good_recovery_rate,
  prevalence        = bad_numeric_input
), expected_error_msg_double)