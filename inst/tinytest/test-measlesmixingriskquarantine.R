# Test just this file: tinytest::run_test_file("inst/tinytest/test-measlesmixingriskquarantine.R")
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

measles_model <- measles::ModelMeaslesMixingRiskQuarantine(
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
# Factory function for creating models with good default values
model_factory <- function(
  n = 9e3,
  prevalence = 1 / 9e3,
  contact_rate = 15,
  transmission_rate = 0.9,
  vax_efficacy = 0.97,
  incubation_period = 10,
  prodromal_period = 3,
  rash_period = 7,
  contact_matrix = cmatrix,
  hospitalization_rate = 0.1,
  hospitalization_period = 10,
  days_undetected = 2,
  quarantine_period_high = 21L,
  quarantine_period_medium = 14L,
  quarantine_period_low = 7L,
  quarantine_willingness = 0.9,
  isolation_willingness = 0.8,
  isolation_period = 10L,
  prop_vaccinated = 0.95,
  detection_rate_quarantine = 0.5,
  contact_tracing_success_rate = 0.8,
  contact_tracing_days_prior = 4L
) {
  ModelMeaslesMixingRiskQuarantine(
    n = n,
    prevalence = prevalence,
    contact_rate = contact_rate,
    transmission_rate = transmission_rate,
    vax_efficacy = vax_efficacy,
    incubation_period = incubation_period,
    prodromal_period = prodromal_period,
    rash_period = rash_period,
    contact_matrix = contact_matrix,
    hospitalization_rate = hospitalization_rate,
    hospitalization_period = hospitalization_period,
    days_undetected = days_undetected,
    quarantine_period_high = quarantine_period_high,
    quarantine_period_medium = quarantine_period_medium,
    quarantine_period_low = quarantine_period_low,
    quarantine_willingness = quarantine_willingness,
    isolation_willingness = isolation_willingness,
    isolation_period = isolation_period,
    prop_vaccinated = prop_vaccinated,
    detection_rate_quarantine = detection_rate_quarantine,
    contact_tracing_success_rate = contact_tracing_success_rate,
    contact_tracing_days_prior = contact_tracing_days_prior
  )
}

# Bad input values
bad_numeric_input <- "not a number"
bad_int_input <- "not an integer"

expected_error_msg_na <- "must not be NA"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

# Test with default values (should work)
expect_silent(model_factory())

# Test each parameter with bad input
expect_error(model_factory(n = bad_int_input), expected_error_msg_int)
expect_error(model_factory(prevalence = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(transmission_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(vax_efficacy = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(incubation_period = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(prodromal_period = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(rash_period = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_matrix = c(1, 0, NA)), expected_error_msg_any_na)
expect_error(model_factory(hospitalization_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(hospitalization_period = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(days_undetected = NA), expected_error_msg_na)
expect_error(model_factory(quarantine_period_high = bad_int_input), expected_error_msg_int)
expect_error(model_factory(quarantine_period_medium = bad_int_input), expected_error_msg_int)
expect_error(model_factory(quarantine_period_low = bad_int_input), expected_error_msg_int)
expect_error(model_factory(quarantine_willingness = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(isolation_willingness = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(isolation_period = bad_int_input), expected_error_msg_int)
expect_error(model_factory(prop_vaccinated = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(detection_rate_quarantine = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_tracing_success_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_tracing_days_prior = bad_int_input), expected_error_msg_int)

# Test bound checks for probability parameters
expect_error(model_factory(vax_efficacy = 1.5), "must be")
expect_error(model_factory(quarantine_willingness = -0.1), "must be")
expect_error(model_factory(detection_rate_quarantine = 1.5), "must be")

