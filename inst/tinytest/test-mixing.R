# Test just this file: tinytest::run_test_file("inst/tinytest/test-mixing.R")

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

flu_model <- ModelSEIRMixing(
  name              = "Flu",
  n                 = N,
  prevalence        = 1 / N,
  contact_rate      = 40,
  transmission_rate = 1.0,
  recovery_rate     = 1 / 10,
  incubation_days   = .009,
  contact_matrix    = cmatrix
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

bad_name <- NA
bad_n <- NA
bad_prevalence <- NA
bad_contact_rate <- NA
bad_transmission_rate <- NA
bad_recovery_rate <- NA
bad_incubation_days <- NA
bad_contact_matrix <- NA

expected_error_msg_str <- "must be a string"
expected_error_msg_int <- "must be an integer"
expected_error_msg_double <- "must be a double"
expected_error_msg_any_na <- "must not contain NA values"

expect_error(test_model <- ModelSEIRMixing(
  name              = bad_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_str)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = bad_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_int)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = bad_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = bad_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = bad_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = bad_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = bad_incubation_days,
  contact_matrix    = good_contact_matrix
), expected_error_msg_double)

expect_error(test_model <- ModelSEIRMixing(
  name              = good_name,
  n                 = good_n,
  prevalence        = good_prevalence,
  contact_rate      = good_contact_rate,
  transmission_rate = good_transmission_rate,
  recovery_rate     = good_recovery_rate,
  incubation_days   = good_incubation_days,
  contact_matrix    = bad_contact_matrix
), expected_error_msg_any_na)