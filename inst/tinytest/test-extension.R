# Test just this file: tinytest::run_test_file("inst/tinytest/test-extension.R")

# C++ extensions (cpp_source_epiworld()): compiling takes minutes, so these
# tests run at home only
if (!at_home() || !requireNamespace("cpp11", quietly = TRUE)) {
  exit_file("C++ extensions are compiled at home only")
}

# An extension reaches the model built in R and adds a global event to it:
# here, one that gives the model's vaccine to every susceptible agent on day 0
cpp_source_epiworld(code = '
class VaccinateAll : public epiworld::GlobalEvent {
public:
  VaccinateAll() { set_name("Vaccinate all"); set_day(0); }
  void operator()(epiworld::Model * m, int) override {
    auto & vaccine = m->get_tool("Vaccine");
    for (auto & agent : m->get_agents())
      if (m->get_stage(agent) == 0)
        agent.add_tool(*m, vaccine);
  }
  std::unique_ptr< epiworld::GlobalEvent > clone_ptr() const override {
    return std::make_unique< VaccinateAll >(*this);
  }
};

[[cpp11::register]]
int vaccinate_all(SEXP model) {
  VaccinateAll event;
  auto & m = epiworldR::as_model(model);
  m.add_globalevent(event);
  return static_cast< int >(m.size());
}
', quiet = TRUE)

new_model <- function() {
  m <- ModelSIR(
    "flu",
    n                 = 1000,
    prevalence        = .01,
    transmission_prob = .5,
    recovery_rate     = .1
  )
  add_layer(m, layer_network("contacts", graph_smallworld(10, .1)))
  add_tool(m, tool("Vaccine", 1, 0, 0, 0), assign_count(0))
  verbose_off(m)
  m
}

susceptible_on_day_30 <- function(m) {
  run(m, ndays = 30, seed = 1)
  totals <- get_hist_total(m)
  sum(totals$counts[totals$date == 30 & totals$stage == "Susceptible"])
}

model <- new_model()
expect_equal(vaccinate_all(model), 1000L)

# The event gives full protection at the end of day 0 to everyone still
# susceptible, so only day 0's infections spread; without it, the epidemic
# runs its course
vaccinated <- susceptible_on_day_30(model)
unvaccinated <- susceptible_on_day_30(new_model())
expect_true(vaccinated > 900)
expect_true(unvaccinated < vaccinated - 300)

# Not a model, and an extension built for another configuration
expect_error(vaccinate_all(1), "not an epiworldR model")
cpp11::cpp_source(code = '
[[cpp11::linking_to("epiworldR")]] void link_();
#define epiworld_double float
#include <epiworldR/extension.h>
[[cpp11::register]] int float_size(SEXP model) {
  return static_cast< int >(epiworldR::as_model(model).size());
}', cxx_std = "CXX17", quiet = TRUE)
expect_error(float_size(model), "Rebuild it against the installed epiworldR")

expect_error(cpp_source_epiworld(), "either")
