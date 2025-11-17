# Test just this file: tinytest::run_test_file("inst/tinytest/test-measles-quarantine-errors.R")

# Factory function for creating models with good default values
model_factory <- function(
  n = 500,
  prevalence = 1,
  contact_rate = 0.1,
  transmission_rate = 0.1,
  vax_efficacy = 0.1,
  incubation_period = 10,
  prodromal_period = 2,
  rash_period = 4,
  days_undetected = 3,
  hospitalization_rate = 0.5,
  hospitalization_period = 5,
  prop_vaccinated = 13 / 15,
  quarantine_period = 15L,
  quarantine_willingness = 1,
  isolation_period = 3L,
  ...
) {
  measles::ModelMeaslesSchool(
    n = n,
    prevalence = prevalence,
    contact_rate = contact_rate,
    transmission_rate = transmission_rate,
    vax_efficacy = vax_efficacy,
    incubation_period = incubation_period,
    prodromal_period = prodromal_period,
    rash_period = rash_period,
    days_undetected = days_undetected,
    hospitalization_rate = hospitalization_rate,
    hospitalization_period = hospitalization_period,
    prop_vaccinated = prop_vaccinated,
    quarantine_period = quarantine_period,
    quarantine_willingness = quarantine_willingness,
    isolation_period = isolation_period,
    ...
  )
}

# Create a MeaslesQuarantine Model --------------------------------------------------------
expect_warning(model_factory(vax_improved_recovery = 0.1), 'removed')
expect_silent(model_factory())
measles_model <- model_factory()

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
# Bad input values
bad_input <- "not a number"
bad_int_input <- "not an integer"

expected_error_msg_na <- "must not be NA"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"

# Test each parameter with bad input
expect_error(model_factory(n = bad_input), expected_error_msg_int)
expect_error(model_factory(prevalence = bad_input), expected_error_msg_int)
expect_error(model_factory(contact_rate = bad_input), expected_error_msg_double)

# Special test for transmission_rate because used in calculation of default contact_rate
expect_error(model_factory(transmission_rate = bad_input), expected_error_msg_double)

expect_error(model_factory(vax_efficacy = bad_input), expected_error_msg_double)
expect_error(model_factory(incubation_period = bad_input), expected_error_msg_double)

# Special test for prodromal_period because used in calculation of default contact_rate
expect_error(model_factory(prodromal_period = bad_input), expected_error_msg_double)

expect_error(model_factory(rash_period = bad_input), expected_error_msg_double)
expect_error(model_factory(days_undetected = bad_input), expected_error_msg_double)
expect_error(model_factory(hospitalization_rate = bad_input), expected_error_msg_double)
expect_error(model_factory(hospitalization_period = bad_input), expected_error_msg_double)
expect_error(model_factory(prop_vaccinated = bad_input), expected_error_msg_double)
expect_error(model_factory(quarantine_period = bad_int_input), expected_error_msg_int)
expect_error(model_factory(quarantine_willingness = bad_input), expected_error_msg_double)
expect_error(model_factory(isolation_period = bad_int_input), expected_error_msg_int)

# Test NA values
expect_error(model_factory(n = NA), expected_error_msg_na)
expect_error(model_factory(prevalence = NA), expected_error_msg_na)
expect_error(model_factory(contact_rate = NA), expected_error_msg_na)
expect_error(model_factory(transmission_rate = NA), expected_error_msg_na)
expect_error(model_factory(vax_efficacy = NA), expected_error_msg_na)
expect_error(model_factory(incubation_period = NA), expected_error_msg_na)
expect_error(model_factory(prodromal_period = NA), expected_error_msg_na)
expect_error(model_factory(rash_period = NA), expected_error_msg_na)
expect_error(model_factory(days_undetected = NA), expected_error_msg_na)
expect_error(model_factory(hospitalization_rate = NA), expected_error_msg_na)
expect_error(model_factory(hospitalization_period = NA), expected_error_msg_na)
expect_error(model_factory(prop_vaccinated = NA), expected_error_msg_na)
expect_error(model_factory(quarantine_period = NA), expected_error_msg_na)
expect_error(model_factory(quarantine_willingness = NA), expected_error_msg_na)
expect_error(model_factory(isolation_period = NA), expected_error_msg_na)