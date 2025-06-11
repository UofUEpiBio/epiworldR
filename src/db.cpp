
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
cpp11::data_frame get_hist_total_cpp(
  SEXP model
) {

  // Making some room
  std::vector< int > date;
  std::vector< std::string > state;
  std::vector< int > counts;

  cpp11::external_pointer<Model<>> ptr(model);
  ptr->get_db().get_hist_total(&date, &state, &counts);

  // Preparing the output
  cpp11::writable::data_frame res({
    "date"_nm   = date,
    "state"_nm  = state,
    "counts"_nm = counts
  });


  return res;
}

[[cpp11::register]]
cpp11::data_frame get_hist_virus_cpp(
  SEXP model
) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector<int> date;
  std::vector<int> id;
  std::vector<std::string> state;
  std::vector<int> counts;

  ptr->get_db().get_hist_virus(
    date, id, state, counts
  );

  // Mapping the id to the name
  std::vector< std::string > viruses;
  for (auto i : ptr->get_viruses())
    viruses.push_back(i->get_name());

  // Mapping using std::transform
  std::vector< std::string > vnames(id.size());
  std::transform(
    id.begin(), id.end(), vnames.begin(),
    [&viruses](int i) { return viruses[i]; }
  );

  return cpp11::writable::data_frame({
    "date"_nm      = date,
    "virus_id"_nm  = id,
    "virus"_nm     = vnames,
    "state"_nm     = state,
    "counts"_nm    = counts,
  });

}

[[cpp11::register]]
cpp11::data_frame get_hist_tool_cpp(
    SEXP model
) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector<int> date;
  std::vector<int> id;
  std::vector<std::string> state;
  std::vector<int> counts;

  ptr->get_db().get_hist_tool(
      date, id, state, counts
  );

  // Same as before, but with tools
  std::vector< std::string > tools;
  for (auto i : ptr->get_viruses())
    tools.push_back(i->get_name());

  std::vector< std::string > tnames(id.size());
  std::transform(
      id.begin(), id.end(), tnames.begin(),
      [&tools](int i) { return tools[i]; }
  );

  return cpp11::writable::data_frame({
    "date"_nm    = date,
    "tool_id"_nm = id,
    "tool"_nm    = tnames,
    "state"_nm   = state,
    "counts"_nm  = counts,
  });

}

[[cpp11::register]]
doubles get_transition_probability_cpp(
    SEXP model
) {

  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::doubles(
    ptr->get_db().get_transition_probability(false)
  );

}

[[cpp11::register]]
cpp11::data_frame get_hist_transition_matrix_cpp(
  SEXP model,
  bool skip_zeros
) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector< std::string > state_from;
  std::vector< std::string > state_to;
  std::vector< int > date;
  std::vector< int > counts;

  ptr->get_db().get_hist_transition_matrix(
    state_from, state_to, date, counts, skip_zeros
  );

  return cpp11::writable::data_frame({
    "state_from"_nm = state_from,
    "state_to"_nm   = state_to,
    "date"_nm       = date,
    "counts"_nm     = counts,
  });

}

[[cpp11::register]]
cpp11::data_frame get_reproductive_number_cpp(
    SEXP model
) {

  // Making some room
  std::vector< int > virus;
  std::vector< int > source;
  std::vector< int > source_exposure_date;
  std::vector< int > counts;

  // Getting the right class
  cpp11::external_pointer<Model<>> ptr(model);
  std::unordered_map< std::vector< int >, int, epiworld::vecHasher<int>> rn =
    ptr->get_db().get_reproductive_number();


  for (const auto & m : rn)
  {
    virus.push_back(m.first[0u]);
    source.push_back(m.first[1u]);
    source_exposure_date.push_back(m.first[2u]);
    counts.push_back(m.second);
  }

  // Same as before: Need to map virus (id) to their names
  std::vector< std::string > viruses;
  for (const auto & i : ptr->get_viruses())
    viruses.push_back(i->get_name());

  std::vector< std::string > vnames(virus.size());
  std::transform(
    virus.begin(), virus.end(), vnames.begin(),
    [&viruses](int i) { return viruses[i]; }
  );

  return cpp11::writable::data_frame({
    "virus_id"_nm           = virus,
    "virus"_nm              = vnames,
    "source"_nm               = source,
    "source_exposure_date"_nm = source_exposure_date,
    "rt"_nm                   = counts
  });

}

[[cpp11::register]]
cpp11::data_frame get_transmissions_cpp(
  SEXP model
) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector<int> date;
  std::vector<int> source;
  std::vector<int> target;
  std::vector<int> virus;
  std::vector<int> source_exposure_date;

  ptr->get_db().get_transmissions(
    date,
    source,
    target,
    virus,
    source_exposure_date
  );

  // Same idea with the names
  std::vector< std::string > viruses;
  for (const auto & i : ptr->get_viruses())
    viruses.push_back(i->get_name());

  std::vector< std::string > vnames(virus.size());
  std::transform(
    virus.begin(), virus.end(), vnames.begin(),
    [&viruses](int i) { return viruses[i]; }
  );

  return cpp11::writable::data_frame({
    "date"_nm                 = date,
    "source"_nm               = source,
    "target"_nm               = target,
    "virus_id"_nm             = virus,
    "virus"_nm                = vnames,
    "source_exposure_date"_nm = source_exposure_date,
  });

}

[[cpp11::register]]
cpp11::data_frame get_generation_time_cpp(
  SEXP model
) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector<int> agent_id;
  std::vector<int> virus_id;
  std::vector<int> date;
  std::vector<int> gentime;

  ptr->get_db().get_generation_time(
    agent_id,
    virus_id,
    date,
    gentime
  );

  // Samething
  std::vector< std::string > viruses;
  for (const auto & i : ptr->get_viruses())
    viruses.push_back(i->get_name());

  std::vector< std::string > vnames(virus_id.size());
  std::transform(
    virus_id.begin(), virus_id.end(), vnames.begin(),
    [&viruses](int i) { return viruses[i]; }
  );

  return cpp11::writable::data_frame({
    "agent"_nm    = agent_id,
    "virus_id"_nm = virus_id,
    "virus"_nm    = vnames,
    "date"_nm     = date,
    "gentime"_nm  = gentime
  });

}

[[cpp11::register]]
cpp11::writable::doubles get_today_total_cpp(SEXP model) {

  cpp11::external_pointer<Model<>> ptr(model);

  std::vector< int > totals;
  std::vector< std::string > names;
  ptr->get_db().get_today_total(&names, &totals);

  cpp11::writable::doubles totals_r(totals.begin(), totals.end());
  cpp11::writable::strings names_r(names.begin(), names.end());

  totals_r.names() = names_r;

  return totals_r;

}
