# Test just this file: tinytest::run_test_file("inst/tinytest/test-diseases.R")

# Several diseases in one model (epiworld NH-5): the flu (created by the
# model) with a second strain, and RSV, a second disease. One infection per
# agent at a time; every output has a disease column.
n <- 3000L
model <- ModelSIR(
  name              = "Flu",
  n                 = n,
  prevalence        = 0.01,
  transmission_prob = 0.3,
  recovery_rate     = 0.2
)
add_layer(model, layer_network("contacts", graph_smallworld(k = 8, p = .05)))
verbose_off(model)

# A strain is the same disease, with its own parameter
add_param(model, "Flu B recovery", 0.4)
flu_b <- strain(get_virus(model, 0), "Flu B")
expect_inherits(flu_b, "epiworld_virus")
use_param(flu_b, "Recovery rate", "Flu B recovery")
add_virus(model, flu_b, assign_count(10))
expect_equal(get_diseases(model), "Flu")

rsv <- virus("RSV") |>
  add_stage("Susceptible") |>
  add_stage("Infected", infectious = TRUE) |>
  add_stage("Recovered") |>
  add_infection("Susceptible", "Infected", rate = 0.2) |>
  add_rate("Infected", "Recovered", rate = 0.1)
add_virus(model, rsv, assign_count(20))
expect_equal(get_diseases(model), c("Flu", "RSV"))

run(model, ndays = 120, seed = 331)

# Each disease's counts add up to the population every day ---------------------
ht <- get_hist_total(model)
expect_equal(unique(ht$disease), c("Flu", "RSV"))
tot <- tapply(ht$counts, list(ht$disease, ht$date), sum)
expect_true(all(tot == n))
expect_equal(unique(get_hist_total(model, disease = "RSV")$disease), "RSV")
expect_error(get_hist_total(model, disease = "Measles"), "no disease \"Measles\"")
expect_error(get_hist_total(model, disease = 1), "disease's name")

# Some agents had both, one after the other
stages_flu <- get_agents_stages(model, disease = "Flu")
stages_rsv <- get_agents_stages(model, disease = "RSV")
expect_true(sum(stages_flu == "Recovered" & stages_rsv == "Recovered") > 0)

# What depends on one disease's stages needs the disease -----------------------
expect_error(get_agents_stages(model), "disease")
expect_error(get_stages(model), "disease")
expect_equal(get_stages(model, disease = "RSV"), c("Susceptible", "Infected", "Recovered"))
expect_error(get_stage(model, get_agents(model)[0]), "disease")
expect_true(get_stage(model, get_agents(model)[0], disease = "Flu") %in% get_stages(model, "Flu"))
expect_error(get_transition_probability(model), "disease")
tp <- get_transition_probability(model, margin = "stage", disease = "RSV")
expect_equal(dim(tp), c(3L, 3L))
expect_equal(dim(get_transition_probability(model, margin = "status")), c(1L, 1L))
expect_error(plot_incidence(model), "several diseases")
expect_silent(inc <- plot_incidence(model, disease = "Flu", plot = FALSE))

# Transmissions and viruses by disease -----------------------------------------
tr <- get_transmissions(model)
expect_equal(sort(unique(tr$disease)), c("Flu", "RSV"))
expect_true(all(tr$virus[tr$disease == "RSV"] == "RSV"))
expect_true(all(tr$virus[tr$disease == "Flu"] %in% c("Flu", "Flu B")))
hv <- get_hist_virus(model, disease = "Flu")
expect_equal(sort(unique(hv$virus)), c("Flu", "Flu B"))

# Transitions by strain add up to the disease's, and each strain recovers at
# its own rate ------------------------------------------------------------------
vt <- get_hist_virus_transitions(model, disease = "Flu")
expect_equal(sort(unique(vt$virus)), c("Flu", "Flu B"))
ht_tr <- get_hist_transitions(model, disease = "Flu")
by_day <- aggregate(counts ~ date + transition, data = vt, FUN = sum)
m <- merge(by_day, ht_tr, by = c("date", "transition"))
expect_equal(nrow(m), nrow(ht_tr))
expect_equal(m$counts.x, m$counts.y)

recovered <- tapply(
  vt$counts[vt$transition == "Infected -> Recovered"],
  vt$virus[vt$transition == "Infected -> Recovered"],
  sum
)
infected_days <- tapply(
  hv$counts[hv$stage == "Infected"],
  hv$virus[hv$stage == "Infected"],
  sum
)
rate <- recovered / infected_days[names(recovered)]
expect_equal(unname(rate["Flu"]), 0.2, tolerance = 0.1)
expect_equal(unname(rate["Flu B"]), 0.4, tolerance = 0.1)

# Every output table has the disease column -------------------------------------
for (out in list(
  get_today_total(model),
  get_hist_tool(model),
  get_hist_transition_matrix(model),
  get_reproductive_number(model),
  get_generation_time(model),
  get_active_cases(model),
  get_outbreak_size(model),
  get_hospitalizations(model)
))
  expect_equal(names(out)[1L], "disease")

# Variants born by mutation are labelled too -----------------------------------
# Mutants are in the database, not among the model's viruses. Mutating needs a
# C++ extension, so this runs at home only.
toolchain <- c(
  "cpp11", "brio", "callr", "cli", "decor", "desc", "glue", "tibble", "vctrs"
)
if (at_home() &&
  all(vapply(toolchain, requireNamespace, logical(1), quietly = TRUE))) {
  cpp_source_epiworld(code = '
  [[cpp11::register]]
  int mutate_every_day(SEXP model) {
    auto & m = epiworldR::as_model(model);
    m.get_virus(0).set_mutation(
      [](epiworld::Agent *, epiworld::Virus &, epiworld::Model *) -> bool {
        return true;
      }
    );
    return 0;
  }
  ', quiet = TRUE)

  m <- ModelSIR("Flu", n = 500, prevalence = 0.05, transmission_prob = 0.4,
    recovery_rate = 0.1)
  add_layer(m, layer_network("contacts", graph_smallworld(k = 8, p = .05)))
  verbose_off(m)
  mutate_every_day(m)
  run(m, ndays = 20, seed = 77)

  hv <- get_hist_virus(m)
  expect_true(length(unique(hv$virus_id)) > 1L)
  rn <- get_reproductive_number(m)
  ht <- get_hist_virus_transitions(m)
  for (out in list(hv, rn, ht)) {
    expect_true(nrow(out) > 0L)
    expect_false(anyNA(out$virus))
    expect_true(all(out$disease == "Flu"))
  }
}
