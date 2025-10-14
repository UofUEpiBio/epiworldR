# Test just this file: tinytest::run_test_file("inst/tinytest/test-measles-quarantine-errors.R")

# Create a MeaslesQuarantine Model --------------------------------------------------------
good_n <- 500
good_prevalence <- 1
good_contact_rate <- 0.1
good_transmission_rate <- 0.1
good_vax_efficacy <- 0.1
good_vax_improved_recovery <- 0.1
good_incubation_period <- 10
good_prodromal_period <- 2
good_rash_period <- 4
good_days_undetected <- 3
good_hospitalization_rate <- 0.5
good_hospitalization_period <- 5
good_prop_vaccinated <- 13 / 15
good_quarantine_period <- 15
good_quarantine_willingness <- 1
good_isolation_period <- 3

expect_silent(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_period = good_hospitalization_period,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_period = good_quarantine_period,
    quarantine_willingness = good_quarantine_willingness,
    isolation_period = good_isolation_period
))

# Check model initialization
expect_inherits(measles_model, "epiworld_measlesquarantine")
expect_inherits(measles_model, "epiworld_model")
expect_length(class(measles_model), 2)

# Check plot
expect_silent(verbose_off(measles_model))
expect_error(plot(measles_model), "model must be run before it can be plotted")
expect_silent(run(measles_model, ndays = 100, seed = 1912))
expect_silent(plot(measles_model))

# Check functions fail with invalid inputs -------------------------------------
bad_input <- "not a number"
bad_prevalence <- NA
bad_contact_rate <- NA
bad_transmission_rate <- NA
bad_vax_efficacy <- NA
bad_vax_improved_recovery <- NA
bad_incubation_period <-  NA
bad_prodromal_period <- NA
bad_rash_period <- NA
bad_days_undetected <- NA
bad_hospitalization_rate <- NA
bad_hospitalization_period <- NA
bad_prop_vaccinated <- NA
bad_quarantine_period <-  NA
bad_quarantine_willingness <- NA
bad_isolation_period <- NA

expected_error_msg_na <- "must not be NA"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

expect_error(measles_model <- ModelMeaslesSchool(
    n = bad_input
), expected_error_msg_int)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prevalence = bad_input
), expected_error_msg_int)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    contact_rate = bad_input
), expected_error_msg_double)

# Special test for transmission_rate because used in calculation of default contact_rate
expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    transmission_rate = bad_input
), "non-numeric argument to binary operator")

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    contact_rate = good_contact_rate,
    transmission_rate = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    vax_efficacy = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    vax_improved_recovery = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    incubation_period = bad_input
), expected_error_msg_double)

# Special test for prodromal_period because used in calculation of default contact_rate
expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prodromal_period = bad_input
), "non-numeric argument to binary operator")

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    contact_rate = good_contact_rate,
    prodromal_period = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    rash_period = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    days_undetected = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    hospitalization_rate = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    hospitalization_period = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prop_vaccinated = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    quarantine_period = bad_input
), expected_error_msg_int)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    quarantine_willingness = bad_input
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    isolation_period = bad_input
), expected_error_msg_int)

# Check handling NA
expect_error(measles_model <- ModelMeaslesSchool(
    n = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prevalence = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    contact_rate = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    transmission_rate = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    vax_efficacy = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    vax_improved_recovery = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    incubation_period = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prodromal_period = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    rash_period = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    days_undetected = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    hospitalization_rate = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    hospitalization_period = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    prop_vaccinated = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    quarantine_period = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    quarantine_willingness = NA
), expected_error_msg_na)

expect_error(measles_model <- ModelMeaslesSchool(
    n = good_n,
    isolation_period = NA
), expected_error_msg_na)