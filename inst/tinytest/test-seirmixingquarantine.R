# Test just this file: tinytest::run_test_file("inst/tinytest/test-seirmixingquarantine.R")

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

flu_model <- ModelSEIRMixingQuarantine(
  name                         = "Flu",
  n                           = N,
  prevalence                  = 1 / N,
  contact_rate                = 40,
  transmission_rate           = 1.0,
  recovery_rate               = 1 / 10,
  incubation_days             = .009,
  contact_matrix              = cmatrix,
  hospitalization_rate        = 0.05,
  hospitalization_period      = 7,
  days_undetected             = 3,
  quarantine_period           = 14,
  quarantine_willingness      = 0.8,
  isolation_willingness       = 0.5,
  isolation_period            = 7,
  contact_tracing_success_rate = 0.7,
  contact_tracing_days_prior  = 3
)

# Adding the entities 
flu_model |>
  add_entity(e1) |>
  add_entity(e2) |>
  add_entity(e3)

run(flu_model, ndays = 100, seed = 1233)
summary(flu_model)

library(data.table)

agents_entities <- lapply(get_entities(flu_model), \(e) {
  entity_get_agents(e)
}) |> rbindlist()

transmissions <- get_transmissions(flu_model) |>
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
  name = "A Virus",
  n = 9e3,
  prevalence = 1 / 9e3,
  contact_rate = 40,
  transmission_rate = 1.0,
  recovery_rate = 1 / 10,
  incubation_days = .009,
  contact_matrix = cmatrix,
  hospitalization_rate = 0.05,
  hospitalization_period = 7,
  days_undetected = 3,
  quarantine_period = 14L,
  quarantine_willingness = 0.8,
  isolation_willingness = 0.5,
  isolation_period = 7L,
  contact_tracing_success_rate = 0.7,
  contact_tracing_days_prior = 3L
) {
  ModelSEIRMixingQuarantine(
    name = name,
    n = n,
    prevalence = prevalence,
    contact_rate = contact_rate,
    transmission_rate = transmission_rate,
    recovery_rate = recovery_rate,
    incubation_days = incubation_days,
    contact_matrix = contact_matrix,
    hospitalization_rate = hospitalization_rate,
    hospitalization_period = hospitalization_period,
    days_undetected = days_undetected,
    quarantine_period = quarantine_period,
    quarantine_willingness = quarantine_willingness,
    isolation_willingness = isolation_willingness,
    isolation_period = isolation_period,
    contact_tracing_success_rate = contact_tracing_success_rate,
    contact_tracing_days_prior = contact_tracing_days_prior
  )
}

# Bad input values
bad_name <- 10
bad_numeric_input <- "not a number"
bad_int_input <- "not an integer"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

# Test with default values (should work)
expect_silent(model_factory())

# Test each parameter with bad input
expect_error(model_factory(name = bad_name), expected_error_msg_str)
expect_error(model_factory(n = bad_numeric_input), expected_error_msg_int)
expect_error(model_factory(prevalence = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(transmission_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(recovery_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(incubation_days = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_matrix = c(1, 0, NA)), expected_error_msg_any_na)
expect_error(model_factory(hospitalization_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(hospitalization_period = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(quarantine_period = bad_int_input), expected_error_msg_int)
expect_error(model_factory(quarantine_willingness = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(isolation_willingness = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(isolation_period = bad_int_input), expected_error_msg_int)
expect_error(model_factory(contact_tracing_success_rate = bad_numeric_input), expected_error_msg_double)
expect_error(model_factory(contact_tracing_days_prior = bad_int_input), expected_error_msg_int)

# Test NA values
expect_error(model_factory(name = NA), expected_error_msg_str)
expect_error(model_factory(days_undetected = NA), expected_error_msg_na)