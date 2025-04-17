# Test just this file: tinytest::run_test_file("inst/tinytest/test-measles-quarantine.R")

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
good_hospitalization_duration <- 5
good_prop_vaccinated <- 13 / 15
good_quarantine_days <- 15
good_quarantine_willingness <- 1

expect_silent(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
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
bad_n <- NA
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
bad_hospitalization_duration <- NA
bad_prop_vaccinated <- NA
bad_quarantine_days <-  NA
bad_quarantine_willingness <- NA

expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = bad_n,
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
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_int)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = bad_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_int)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = bad_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = bad_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = bad_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = bad_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = bad_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = bad_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = bad_rash_period,
    days_undetected = good_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
    n = good_n,
    prevalence = good_prevalence,
    contact_rate = good_contact_rate,
    transmission_rate = good_transmission_rate,
    vax_efficacy = good_vax_efficacy,
    vax_improved_recovery = good_vax_improved_recovery,
    incubation_period = good_incubation_period,
    prodromal_period = good_prodromal_period,
    rash_period = good_rash_period,
    days_undetected = bad_days_undetected,
    hospitalization_rate = good_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_rate = bad_hospitalization_rate,
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_duration = bad_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = bad_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = bad_quarantine_days,
    quarantine_willingness = good_quarantine_willingness
), expected_error_msg_double)

expect_error(measles_model <- ModelMeaslesQuarantine(
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
    hospitalization_duration = good_hospitalization_duration,
    prop_vaccinated = good_prop_vaccinated,
    quarantine_days = good_quarantine_days,
    quarantine_willingness = bad_quarantine_willingness
), expected_error_msg_double)