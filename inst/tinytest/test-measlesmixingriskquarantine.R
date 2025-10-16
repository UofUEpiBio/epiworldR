# Test just this file: tinytest::run_test_file("inst/tinytest/test-measlesmixingriskquarantine.R")

library(epiworldR)

e1 <- entity("Population 1", 3e3, FALSE)
e2 <- entity("Population 2", 3e3, FALSE)
e3 <- entity("Population 3", 3e3, FALSE)

# Row-stochastic matrix (rowsums 1)
cmatrix <- c(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1)
) |> as.double() |> matrix(byrow = TRUE, nrow = 3)

N <- 9e3

measles_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = N,
  prevalence                 = 1 / N,
  contact_rate               = 15,
  transmission_rate          = 0.9,
  vax_efficacy               = 0.97,
  incubation_period          = 10,
  prodromal_period           = 3,
  rash_period                = 7,
  contact_matrix             = cmatrix,
  hospitalization_rate       = 0.1,
  hospitalization_period     = 10,
  days_undetected            = 2,
  quarantine_period_high     = 21,
  quarantine_period_medium   = 14,
  quarantine_period_low      = 7,
  quarantine_willingness     = 0.9,
  isolation_willingness      = 0.8,
  isolation_period           = 10,
  prop_vaccinated            = 0.95,
  detection_rate_quarantine  = 0.5,
  contact_tracing_success_rate = 0.8,
  contact_tracing_days_prior = 4
)

# Adding the entities 
measles_model |>
  add_entity(e1) |>
  add_entity(e2) |>
  add_entity(e3)

run(measles_model, ndays = 100, seed = 1233)
summary(measles_model)

library(data.table)

agents_entities <- lapply(get_entities(measles_model), \(e) {
  entity_get_agents(e)
}) |> rbindlist()

transmissions <- get_transmissions(measles_model) |>
  data.table()

# We only need the date and the source
transmissions <- subset(
  transmissions,
  select = c("date", "target")
)

# Attaching the entity to the source
transmissions <- merge(
  transmissions,
  agents_entities,
  by.x = "target", by.y = "agent"
)

# Aggregating by date x entity (counts)
transmissions <- transmissions[, .N, by = .(date, entity)]

transmissions[, sum(N), by = .(entity)]

transmissions[, table(entity)]

# Check functions fail with invalid inputs -------------------------------------
good_n <- 9e3
good_prevalence <- 1 / good_n
good_contact_rate <- 15
good_transmission_rate <- 0.9
good_vax_efficacy <- 0.97
good_incubation_period <- 10
good_prodromal_period <- 3
good_rash_period <- 7
good_contact_matrix <- cmatrix
good_hospitalization_rate <- 0.1
good_hospitalization_period <- 10
good_days_undetected <- 2
good_quarantine_period_high <- 21L
good_quarantine_period_medium <- 14L
good_quarantine_period_low <- 7L
good_quarantine_willingness <- 0.9
good_isolation_willingness <- 0.8
good_isolation_period <- 10L
good_prop_vaccinated <- 0.95
good_detection_rate_quarantine <- 0.5
good_contact_tracing_success_rate <- 0.8
good_contact_tracing_days_prior <- 4L

bad_numeric_input <- "not a number"
bad_int_input <- "not an integer"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = bad_numeric_input,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_int)

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = bad_numeric_input,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_double)

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = c(1, 0, NA),
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_any_na)

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = bad_numeric_input,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_double)

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = bad_int_input,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_int)

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = NA,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), expected_error_msg_na)

# Test bound checks for probability parameters
expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = 1.5, # Invalid: > 1
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), "must be")

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = -0.1, # Invalid: < 0
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = good_detection_rate_quarantine,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), "must be")

expect_error(test_model <- ModelMeaslesMixingRiskQuarantine(
  n                          = good_n,
  prevalence                 = good_prevalence,
  contact_rate               = good_contact_rate,
  transmission_rate          = good_transmission_rate,
  vax_efficacy               = good_vax_efficacy,
  incubation_period          = good_incubation_period,
  prodromal_period           = good_prodromal_period,
  rash_period                = good_rash_period,
  contact_matrix             = good_contact_matrix,
  hospitalization_rate       = good_hospitalization_rate,
  hospitalization_period     = good_hospitalization_period,
  days_undetected            = good_days_undetected,
  quarantine_period_high     = good_quarantine_period_high,
  quarantine_period_medium   = good_quarantine_period_medium,
  quarantine_period_low      = good_quarantine_period_low,
  quarantine_willingness     = good_quarantine_willingness,
  isolation_willingness      = good_isolation_willingness,
  isolation_period           = good_isolation_period,
  prop_vaccinated            = good_prop_vaccinated,
  detection_rate_quarantine  = 1.5, # Invalid: > 1
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior = good_contact_tracing_days_prior
), "must be")
