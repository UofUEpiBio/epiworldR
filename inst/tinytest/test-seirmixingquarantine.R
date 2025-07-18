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
good_name <- "A Virus"
good_n <- 9e3
good_prevalence <- 1 / good_n
good_contact_rate <- 40
good_transmission_rate <- 1.0
good_recovery_rate <- 1 / 10
good_incubation_days <- .009
good_contact_matrix <- cmatrix
good_hospitalization_rate <- 0.05
good_hospitalization_period <- 7
good_days_undetected <- 3
good_quarantine_period <- 14L
good_quarantine_willingness <- 0.8
good_isolation_willingness <- 0.5
good_isolation_period <- 7L
good_contact_tracing_success_rate <- 0.7
good_contact_tracing_days_prior <- 3L

bad_name <- 10
bad_numeric_input <- "not a number"
bad_int_input <- "not an integer"

expected_error_msg_na <- "must not be NA"
expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = bad_name,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_str)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = bad_numeric_input,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_int)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = good_n,
  prevalence                  = bad_numeric_input,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = c(1, 0, NA),
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_any_na)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = bad_numeric_input,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = bad_int_input,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_int)

# Check NA values
expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = NA,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = good_days_undetected,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_str)

expect_error(test_model <- ModelSEIRMixingQuarantine(
  name                         = good_name,
  n                           = good_n,
  prevalence                  = good_prevalence,
  contact_rate                = good_contact_rate,
  transmission_rate           = good_transmission_rate,
  recovery_rate               = good_recovery_rate,
  incubation_days             = good_incubation_days,
  contact_matrix              = good_contact_matrix,
  hospitalization_rate        = good_hospitalization_rate,
  hospitalization_period      = good_hospitalization_period,
  days_undetected             = NA,
  quarantine_period           = good_quarantine_period,
  quarantine_willingness      = good_quarantine_willingness,
  isolation_willingness       = good_isolation_willingness,
  isolation_period            = good_isolation_period,
  contact_tracing_success_rate = good_contact_tracing_success_rate,
  contact_tracing_days_prior  = good_contact_tracing_days_prior
), expected_error_msg_na)