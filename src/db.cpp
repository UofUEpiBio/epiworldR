
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
  std::vector< std::string > status;
  std::vector< int > counts;
  
  cpp11::external_pointer<Model<>> ptr(model);
  ptr->get_db().get_hist_total(&dates, &status, &counts);
  
  
  // Preparing the output
  cpp11::writable::data_frame res({
    "dates"_nm = dates,
    "status"_nm = status,
    "counts"_nm = counts
  });
  
  
  return res;
}

[[cpp11::register]]
doubles get_transition_probability_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
  
}

[[cpp11::register]]
cpp11::strings get_status_cpp(
    SEXP model
) {
  
  cpp11::external_pointer<Model<>> ptr(model);
  return cpp11::writable::strings(ptr->get_status());
  
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

