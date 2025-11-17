# Test just this file: tinytest::run_test_file("inst/tinytest/test-measlesmixingriskquarantine.R")
# Row-stochastic matrix (rowsums 1)
self <- .83
others <- (1 - self) / 2
cmatrix <- c(
  c(self, others, others),
  c(others, self, others),
  c(others, others, self)
) |> as.double() |> matrix(byrow = TRUE, nrow = 3)

N <- 600

r0 <- 15
c_rate <- 20
i_rate <- r0 / c_rate * (1/3)

model_factory <- function(durations, nsims = 100) {

  e1 <- entity("Population 1", 200, FALSE)
  e2 <- entity("Population 2", 200, FALSE)
  e3 <- entity("Population 3", 200, FALSE)

  measles_model <- ModelMeaslesMixingRiskQuarantine(
    n                          = N,
    prevalence                 = 1 / N,
    contact_rate               = c_rate,
    transmission_rate          = i_rate,
    vax_efficacy               = 0.97,
    incubation_period          = 10,
    prodromal_period           = 3,
    rash_period                = 7,
    contact_matrix             = cmatrix,
    hospitalization_rate       = 0.1,
    hospitalization_period     = 10,
    days_undetected            = 2,
    quarantine_period_high     = durations[1],
    quarantine_period_medium   = durations[2],
    quarantine_period_low      = durations[3],
    quarantine_willingness     = 0.9,
    isolation_willingness      = 0.8,
    isolation_period           = 10,
    prop_vaccinated            = 0.3,
    detection_rate_quarantine  = 0.5,
    contact_tracing_success_rate = 0.8,
    contact_tracing_days_prior = 4
  )

  # Adding the entities 
  measles_model |>
    add_entity(e1) |>
    add_entity(e2) |>
    add_entity(e3)

  run_multiple(
    measles_model, ndays = 60, nsims = nsims, seed = 221,
    saver = make_saver("total_hist"),
    nthreads = 2
    )

  res <- run_multiple_get_results(
    measles_model, freader = data.table::fread,
    nthreads = 1
    )$total_hist

  res[
    (date == get_ndays(measles_model)) &
    !(state %in% c("Susceptible", "Quarantined Susceptible"))
    ][, sum(counts), by = "sim_num"]$V1
}

nsims <- 50
ans_base <- model_factory(c(21L, 21L, 21L), nsims = nsims) |> mean()
ans_none <- model_factory(c(0L, 0L, 0L), nsims = nsims) |> mean()
ans_high <- model_factory(c(21L, 0L, 0L), nsims = nsims) |> mean()
ans_mid  <- model_factory(c(0L, 21L, 0L), nsims = nsims) |> mean()
ans_low  <- model_factory(c(0L, 0L, 21L), nsims = nsims) |> mean()

# Basic expectations. No protection should be the worst
expect_true(ans_none > ans_base)
expect_true(ans_none > ans_high)
expect_true(ans_none > ans_mid)
expect_true(ans_none > ans_low)

# The best protection should be the full quarantine
expect_true(ans_base < ans_high)
expect_true(ans_base < ans_mid)
expect_true(ans_base < ans_low)


