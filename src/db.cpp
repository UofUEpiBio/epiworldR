
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
cpp11::data_frame get_hist_total_cpp(
  SEXP model
) {
  
  // Making some room
  std::vector< int > dates;
  std::vector< std::string > state;
  std::vector< int > counts;
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->get_db().get_hist_total(&dates, &state, &counts);
  
  
  // Preparing the output
  cpp11::writable::data_frame res({
    "dates"_nm = dates,
    "state"_nm = state,
    "counts"_nm = counts
  });
  
  
  return res;
}

[[cpp11::register]]
cpp11::data_frame get_hist_variant_cpp(
  SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  
  std::vector<int> date;
  std::vector<int> id;
  std::vector<std::string> state;
  std::vector<int> counts;
  
  ptr->get_db().get_hist_variant(
    date, id, state, counts
  );
  
  return cpp11::writable::data_frame({
    "date"_nm   = date, 
    "id"_nm     = id,
    "state"_nm  = state,
    "counts"_nm = counts,
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
  
  return cpp11::writable::data_frame({
    "date"_nm   = date, 
      "id"_nm     = id,
      "state"_nm  = state,
      "counts"_nm = counts,
  });
  
}

[[cpp11::register]]
doubles get_transition_probability_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
  
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
  std::vector< int > variant;
  std::vector< int > source;
  std::vector< int > source_exposure_dates;
  std::vector< int > counts;
  
  // Getting the right class
  cpp11::external_pointer<Model<>> ptr(model);
  std::unordered_map< std::vector< int >, int, epiworld::vecHasher<int>> rn =
    ptr->get_db().reproductive_number();
  

  for (const auto & m : rn) 
  {
    variant.push_back(m.first[0u]);
    source.push_back(m.first[1u]);
    source_exposure_dates.push_back(m.first[2u]);
    counts.push_back(m.second);
  }
  
  return cpp11::writable::data_frame({
    "variant"_nm                = variant,
      "source"_nm                = source,
      "source_exposure_dates"_nm = source_exposure_dates,
      "counts"_nm                = counts
  });
  
}

