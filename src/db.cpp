
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

// Every output takes a disease (its index, see get_diseases(); -1: every
// disease, one after the other) and returns a `disease` column first.

// The names of the viruses (variants) with these ids: from the database, which
// knows the variants born by mutation, not only the model's viruses
static std::vector< std::string > virus_names(
  Model & model,
  const std::vector< int > & id
) {

  const auto & db = model.get_db();
  std::vector< std::string > out;
  out.reserve(id.size());
  for (int i : id)
    out.push_back(db.get_virus_name(i));
  return out;

}

[[cpp11::register]]
cpp11::data_frame get_hist_total_cpp(SEXP model, int disease) {

  std::vector< int > date, counts;
  std::vector< std::string > status, stage, dnames;

  cpp11::external_pointer<Model> ptr(model);
  ptr->get_db().get_hist_total(&date, &status, &stage, &counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm = dnames,
    "date"_nm    = date,
    "status"_nm  = status,
    "stage"_nm   = stage,
    "counts"_nm  = counts
  });

}

[[cpp11::register]]
cpp11::data_frame get_hist_virus_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, id, counts;
  std::vector< std::string > status, stage, dnames;
  ptr->get_db().get_hist_virus(date, id, status, stage, counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm  = dnames,
    "date"_nm     = date,
    "virus_id"_nm = id,
    "virus"_nm    = virus_names(*ptr, id),
    "status"_nm   = status,
    "stage"_nm    = stage,
    "counts"_nm   = counts,
  });

}

[[cpp11::register]]
cpp11::data_frame get_hist_tool_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, id, counts;
  std::vector< std::string > status, stage, dnames;
  ptr->get_db().get_hist_tool(date, id, status, stage, counts, disease, &dnames);

  std::vector< std::string > tools;
  for (const auto & t : ptr->get_tools())
    tools.push_back(t->get_name());

  std::vector< std::string > tnames(id.size());
  std::transform(
    id.begin(), id.end(), tnames.begin(),
    [&tools](int i) { return tools[static_cast< size_t >(i)]; }
  );

  return cpp11::writable::data_frame({
    "disease"_nm = dnames,
    "date"_nm    = date,
    "tool_id"_nm = id,
    "tool"_nm    = tnames,
    "status"_nm  = status,
    "stage"_nm   = stage,
    "counts"_nm  = counts,
  });

}

// The transition matrix over a margin ("pair", "stage" or "status"): the
// values (column-major, from x to) and the labels
[[cpp11::register]]
cpp11::list get_transition_probability_cpp(
  SEXP model,
  std::string margin,
  int disease,
  bool normalize
) {

  cpp11::external_pointer<Model> ptr(model);

  epiworld::Margin m;
  if (margin == "pair")
    m = epiworld::Margin::pair;
  else if (margin == "stage")
    m = epiworld::Margin::stage;
  else if (margin == "status")
    m = epiworld::Margin::status;
  else
    cpp11::stop(
      "-margin- must be \"pair\", \"stage\" or \"status\" (got \"%s\").",
      margin.c_str()
    );

  auto tm = ptr->get_db().get_transition_probability(m, disease, normalize);

  return cpp11::writable::list({
    "values"_nm = cpp11::writable::doubles(tm.values),
    "status"_nm = cpp11::writable::strings(tm.status),
    "stage"_nm  = cpp11::writable::strings(tm.stage)
  });

}

[[cpp11::register]]
cpp11::data_frame get_hist_transition_matrix_cpp(
  SEXP model,
  bool skip_zeros,
  int disease
) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< std::string > status_from, stage_from, status_to, stage_to, dnames;
  std::vector< int > date, counts;

  ptr->get_db().get_hist_transition_matrix(
    status_from, stage_from, status_to, stage_to, date, counts, skip_zeros,
    disease, &dnames
  );

  return cpp11::writable::data_frame({
    "disease"_nm     = dnames,
    "status_from"_nm = status_from,
    "stage_from"_nm  = stage_from,
    "status_to"_nm   = status_to,
    "stage_to"_nm    = stage_to,
    "date"_nm        = date,
    "counts"_nm      = counts,
  });

}

// How many times each declared transition fired, by day
[[cpp11::register]]
cpp11::data_frame get_hist_transitions_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, counts;
  std::vector< std::string > transition, dnames;
  ptr->get_db().get_hist_transitions(date, transition, counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm    = dnames,
    "date"_nm       = date,
    "transition"_nm = transition,
    "counts"_nm     = counts
  });

}

// The same, by variant (the stage transitions an agent with the virus takes)
[[cpp11::register]]
cpp11::data_frame get_hist_virus_transitions_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, id, counts;
  std::vector< std::string > transition, dnames;
  ptr->get_db().get_hist_virus_transitions(
    date, id, transition, counts, disease, &dnames
  );

  return cpp11::writable::data_frame({
    "disease"_nm    = dnames,
    "date"_nm       = date,
    "virus_id"_nm   = id,
    "virus"_nm      = virus_names(*ptr, id),
    "transition"_nm = transition,
    "counts"_nm     = counts
  });

}

[[cpp11::register]]
cpp11::data_frame get_reproductive_number_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);
  auto rn = ptr->get_db().get_reproductive_number(disease);

  std::vector< int > virus, source, source_exposure_date, counts;
  for (const auto & m : rn)
  {
    virus.push_back(m.first[0u]);
    source.push_back(m.first[1u]);
    source_exposure_date.push_back(m.first[2u]);
    counts.push_back(m.second);
  }

  // The disease of each row: its virus's
  const auto diseases = ptr->get_diseases();
  std::vector< std::string > dnames(virus.size());
  for (size_t i = 0u; i < virus.size(); ++i)
    dnames[i] = diseases[static_cast< size_t >(
      ptr->get_db().get_virus_disease(virus[i])
    )];

  return cpp11::writable::data_frame({
    "disease"_nm              = dnames,
    "virus_id"_nm             = virus,
    "virus"_nm                = virus_names(*ptr, virus),
    "source"_nm               = source,
    "source_exposure_date"_nm = source_exposure_date,
    "rt"_nm                   = counts
  });

}

[[cpp11::register]]
cpp11::data_frame get_transmissions_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, source, target, virus, source_exposure_date, layer;
  std::vector< std::string > dnames;
  ptr->get_db().get_transmissions(
    date, source, target, virus, source_exposure_date, layer, disease, &dnames
  );

  return cpp11::writable::data_frame({
    "disease"_nm              = dnames,
    "date"_nm                 = date,
    "source"_nm               = source,
    "target"_nm               = target,
    "virus_id"_nm             = virus,
    "virus"_nm                = virus_names(*ptr, virus),
    "source_exposure_date"_nm = source_exposure_date,
    "layer"_nm                = layer
  });

}

[[cpp11::register]]
cpp11::data_frame get_generation_time_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > agent_id, virus_id, date, gentime;
  std::vector< std::string > dnames;
  ptr->get_db().get_generation_time(
    agent_id, virus_id, date, gentime, disease, &dnames
  );

  return cpp11::writable::data_frame({
    "disease"_nm  = dnames,
    "agent"_nm    = agent_id,
    "virus_id"_nm = virus_id,
    "virus"_nm    = virus_names(*ptr, virus_id),
    "date"_nm     = date,
    "gentime"_nm  = gentime
  });

}

[[cpp11::register]]
cpp11::writable::data_frame get_today_total_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > counts;
  std::vector< std::string > status, stage, dnames;
  ptr->get_db().get_today_total(&status, &stage, &counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm = dnames,
    "status"_nm  = status,
    "stage"_nm   = stage,
    "counts"_nm  = counts
  });

}

[[cpp11::register]]
cpp11::writable::data_frame get_active_cases_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, virus, counts;
  std::vector< std::string > dnames;
  ptr->get_db().get_active_cases(date, virus, counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm      = dnames,
    "date"_nm         = date,
    "virus_id"_nm     = virus,
    "virus"_nm        = virus_names(*ptr, virus),
    "active_cases"_nm = counts,
  });

}

[[cpp11::register]]
cpp11::writable::data_frame get_outbreak_size_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, virus, counts;
  std::vector< std::string > dnames;
  ptr->get_db().get_outbreak_size(date, virus, counts, disease, &dnames);

  return cpp11::writable::data_frame({
    "disease"_nm       = dnames,
    "date"_nm          = date,
    "virus_id"_nm      = virus,
    "virus"_nm         = virus_names(*ptr, virus),
    "outbreak_size"_nm = counts,
  });

}

[[cpp11::register]]
cpp11::writable::data_frame get_hospitalizations_cpp(SEXP model, int disease) {

  cpp11::external_pointer<Model> ptr(model);

  std::vector< int > date, virus_id, tool_id, count;
  std::vector< double > weight;
  std::vector< std::string > dnames;
  ptr->get_db().get_hospitalizations(
    date, virus_id, tool_id, count, weight, disease, &dnames
  );

  return cpp11::writable::data_frame({
    "disease"_nm  = dnames,
    "date"_nm     = date,
    "virus_id"_nm = virus_id,
    "tool_id"_nm  = tool_id,
    "count"_nm    = count,
    "weight"_nm   = weight
  });

}

// The model's diseases, by index
[[cpp11::register]]
cpp11::writable::strings get_diseases_cpp(SEXP model) {
  cpp11::external_pointer<Model> ptr(model);
  return cpp11::writable::strings(ptr->get_diseases());
}
