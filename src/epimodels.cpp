
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/matrix.hpp"
#include "epiworld-common.h"

using namespace cpp11;
using namespace epiworld;

// Model definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/models

// The registry: every registered model is built by make_model() ------------

// The core models (the measles package registers its own)
static const ModelRegistry & registry() {

  static const ModelRegistry reg = [] {
    ModelRegistry r;
    epimodels::register_models(r);
    return r;
  }();

  return reg;

}

static cpp11::writable::data_frame params_df(const ModelInfo & info) {

  cpp11::writable::strings arg, name, description;
  cpp11::writable::doubles value, min, max;
  cpp11::writable::logicals integer;
  for (const auto & p : info.params)
  {
    arg.push_back(p.arg);
    name.push_back(p.name);
    value.push_back(p.value);
    min.push_back(p.min);
    max.push_back(p.max);
    integer.push_back(p.integer ? TRUE : FALSE);
    description.push_back(p.description);
  }

  return cpp11::writable::data_frame({
    "arg"_nm         = arg,
    "name"_nm        = name,
    "value"_nm       = value,
    "min"_nm         = min,
    "max"_nm         = max,
    "integer"_nm     = integer,
    "description"_nm = description
  });

}

[[cpp11::register]]
cpp11::writable::data_frame list_models_cpp() {

  cpp11::writable::strings id, label, description, virus_name;
  for (const auto & info : registry().list_models())
  {
    id.push_back(info.id);
    label.push_back(info.label);
    description.push_back(info.description);
    virus_name.push_back(info.virus_name);
  }

  return cpp11::writable::data_frame({
    "id"_nm          = id,
    "label"_nm       = label,
    "description"_nm = description,
    "virus_name"_nm  = virus_name
  });

}

[[cpp11::register]]
cpp11::writable::list get_model_info_cpp(std::string id) {

  const ModelInfo & info = registry().get_model_info(id);
  return cpp11::writable::list({
    "id"_nm          = info.id,
    "label"_nm       = info.label,
    "description"_nm = info.description,
    "virus_name"_nm  = info.virus_name,
    "params"_nm      = params_df(info)
  });

}

[[cpp11::register]]
SEXP make_model_cpp(std::string id, cpp11::doubles values, std::string virus_name) {

  std::map< std::string, double > v;
  if (values.size() > 0)
  {
    cpp11::strings names(values.names());
    for (R_xlen_t i = 0; i < values.size(); ++i)
      v[std::string(names[i])] = values[i];
  }

  return cpp11::external_pointer<Model>(
    registry().make_model(id, v, virus_name).release()
  );

}

[[cpp11::register]]
SEXP ModelSIRLogit_cpp(
  std::string vname,
  unsigned int n,
  SEXP data,
  int ncols,
  std::vector< double > coefs_infect,
  std::vector< double > coefs_recover,
  std::vector< int > coef_infect_cols,
  std::vector< int > coef_recover_cols,
  double prob_infection,
  double recovery_rate,
  double prevalence
) {

  std::vector< size_t > cinfect;
  std::vector< size_t > crecover;

  for (auto i : coef_infect_cols)
    cinfect.push_back(static_cast<size_t>(i));

  for (auto i : coef_recover_cols)
    crecover.push_back(static_cast<size_t>(i));

  cpp11::external_pointer<epiworld::epimodels::ModelSIRLogit> ptr(
    new epiworld::epimodels::ModelSIRLogit(
        vname,
        n,
        REAL(data),
        static_cast<size_t>(ncols),
        coefs_infect,
        coefs_recover,
        cinfect,
        crecover,
        prob_infection,
        recovery_rate,
        prevalence
    )
  );

  return ptr;

}

[[cpp11::register]]
SEXP ModelDiffNet_cpp(
  std::string name,
  unsigned int n,
  double prevalence,
  double prob_adopt,
  bool normalize_exposure,
  SEXP data,
  int data_ncols,
  std::vector< int > data_cols,
  std::vector<double> params
) {

  // Maps data_cols to size_t
  std::vector< size_t > data_cols_s;
  for (auto i : data_cols)
    data_cols_s.push_back(static_cast<size_t>(i));

  cpp11::external_pointer<epiworld::epimodels::ModelDiffNet> ptr(
    new epiworld::epimodels::ModelDiffNet(
      name,
      n,
      prevalence,
      prob_adopt,
      normalize_exposure,
      &(REAL(data)[0u]),
      data_ncols,
      data_cols_s,
      params
    )
  );

  return ptr;

}
