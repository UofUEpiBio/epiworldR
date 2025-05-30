# Test just this file: tinytest::run_test_file("inst/tinytest/test-surv.R")

# Create a SURV Model ----------------------------------------------------------
good_name                   <- "A Virus"
good_prevalence             <- 20
good_efficacy_vax           <- 0.6
good_latent_period          <- 4
good_infect_period          <- 5
good_prob_symptoms          <- 0.5
good_prop_vaccinated        <- 0.7
good_prop_vax_redux_transm  <- 0.8
good_prop_vax_redux_infect  <- 0.95
good_surveillance_prob      <- 0.1
good_transmission_rate      <- 0.2
good_prob_death             <- 0.001
good_prob_noreinfect        <- 0.5

expect_silent(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
))

# Check model initialization ---------------------------------------------------
expect_inherits(test_model, "epiworld_surv")
expect_inherits(test_model, "epiworld_model")
expect_length(class(test_model), 2)

# Check functions fail with invalid inputs -------------------------------------
bad_name                <- 10
bad_numeric_input       <- "not a number"

expected_error_msg_str    <- "must be a string"
expected_error_msg_double <- "must be a double"

expect_error(test_model <- ModelSURV(
    name                    = bad_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_str)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = bad_numeric_input,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = bad_numeric_input,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = bad_numeric_input,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = bad_numeric_input,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = bad_numeric_input,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = bad_numeric_input,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = bad_numeric_input,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = bad_numeric_input,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = bad_numeric_input,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = bad_numeric_input,
    prob_death              = good_prob_death,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = bad_numeric_input,
    prob_noreinfect         = good_prob_noreinfect
), expected_error_msg_double)

expect_error(test_model <- ModelSURV(
    name                    = good_name,
    prevalence              = good_prevalence,
    efficacy_vax            = good_efficacy_vax,
    latent_period           = good_latent_period,
    infect_period           = good_infect_period,
    prob_symptoms           = good_prob_symptoms,
    prop_vaccinated         = good_prop_vaccinated,
    prop_vax_redux_transm   = good_prop_vax_redux_transm,
    prop_vax_redux_infect   = good_prop_vax_redux_infect,
    surveillance_prob       = good_surveillance_prob,
    transmission_rate       = good_transmission_rate,
    prob_death              = good_prob_death,
    prob_noreinfect         = bad_numeric_input
), expected_error_msg_double)


