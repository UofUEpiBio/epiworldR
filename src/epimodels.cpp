
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/matrix.hpp"
#include "epiworld-common.h"

using namespace cpp11;

// Measles Model definitions
// Based on epiworld library: https://github.com/UofUEpiBio/epiworld

[[cpp11::register]]
SEXP ModelMeaslesSchool_cpp(
  unsigned int n,
  unsigned int prevalence,
  double contact_rate,
  double transmission_rate,
  double vax_efficacy,
  double vax_reduction_recovery_rate,
  double incubation_period,
  double prodromal_period,
  double rash_period,
  double days_undetected,
  double hospitalization_rate,
  double hospitalization_period,
  double prop_vaccinated,
  int quarantine_period,
  double quarantine_willingness,
  int isolation_period
) {

  // Creating a pointer to a ModelMeaslesSchool model
  cpp11::external_pointer<epiworld::epimodels::ModelMeaslesSchool<>> ptr(
      new epiworld::epimodels::ModelMeaslesSchool<>(
          n,
          prevalence,
          contact_rate,
          transmission_rate,
          vax_efficacy,
          vax_reduction_recovery_rate,
          incubation_period,
          prodromal_period,
          rash_period,
          days_undetected,
          hospitalization_rate,
          hospitalization_period,
          prop_vaccinated,
          quarantine_period,
          quarantine_willingness,
          isolation_period
      )
  );

  return ptr;

}

[[cpp11::register]]
SEXP ModelMeaslesMixing_cpp(
    unsigned int n,
    double prevalence,
    double contact_rate,
    double transmission_rate,
    double vax_efficacy,
    double vax_reduction_recovery_rate,
    double incubation_period,
    double prodromal_period,
    double rash_period,
    std::vector< double > contact_matrix,
    double hospitalization_rate,
    double hospitalization_period,
    // Policy parameters
    double days_undetected,
    int quarantine_period,
    double quarantine_willingness,
    double isolation_willingness,
    int isolation_period,
    double prop_vaccinated,
    double contact_tracing_success_rate = 1.0,
    unsigned int contact_tracing_days_prior = 4u
) {

  // Creating a pointer to a ModelMeaslesMixing model
  cpp11::external_pointer<epiworld::epimodels::ModelMeaslesMixing<>> ptr(
      new epiworld::epimodels::ModelMeaslesMixing<>(
          n,
          prevalence,
          contact_rate,
          transmission_rate,
          vax_efficacy,
          vax_reduction_recovery_rate,
          incubation_period,
          prodromal_period,
          rash_period,
          contact_matrix,
          hospitalization_rate,
          hospitalization_period,
          days_undetected,
          quarantine_period,
          quarantine_willingness,
          isolation_willingness,
          isolation_period,
          prop_vaccinated,
          contact_tracing_success_rate,
          contact_tracing_days_prior
      )
  );

  return ptr;

}

[[cpp11::register]]
SEXP ModelMeaslesMixingRiskQuarantine_cpp(
    unsigned int n,
    double prevalence,
    double contact_rate,
    double transmission_rate,
    double vax_efficacy,
    double incubation_period,
    double prodromal_period,
    double rash_period,
    std::vector< double > contact_matrix,
    double hospitalization_rate,
    double hospitalization_period,
    // Policy parameters
    double days_undetected,
    int quarantine_period_high,
    int quarantine_period_medium,
    int quarantine_period_low,
    double quarantine_willingness,
    double isolation_willingness,
    int isolation_period,
    double prop_vaccinated,
    double detection_rate_quarantine,
    double contact_tracing_success_rate = 1.0,
    unsigned int contact_tracing_days_prior = 4u
) {

  // Creating a pointer to a ModelMeaslesMixingRiskQuarantine model
  cpp11::external_pointer<epiworld::epimodels::ModelMeaslesMixingRiskQuarantine<>> ptr(
      new epiworld::epimodels::ModelMeaslesMixingRiskQuarantine<>(
          n,
          prevalence,
          contact_rate,
          transmission_rate,
          vax_efficacy,
          incubation_period,
          prodromal_period,
          rash_period,
          contact_matrix,
          hospitalization_rate,
          hospitalization_period,
          days_undetected,
          quarantine_period_high,
          quarantine_period_medium,
          quarantine_period_low,
          quarantine_willingness,
          isolation_willingness,
          isolation_period,
          prop_vaccinated,
          detection_rate_quarantine,
          contact_tracing_success_rate,
          contact_tracing_days_prior
      )
  );

  return ptr;

}
