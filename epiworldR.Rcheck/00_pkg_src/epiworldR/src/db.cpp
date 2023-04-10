
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"

using namespace cpp11;

#define WrapModel(model, name) \
  cpp11::external_pointer<epiworld::epimodels::model<>> (name) 

[[cpp11::register]]
int queuing_on_cpp(
  SEXP model,
  std::string model_class
  ) {

  // Getting the right class
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    ptr->queuing_on();
  
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    ptr->queuing_on();
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    ptr->queuing_on();    

  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    ptr->queuing_on();

  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    ptr->queuing_on();
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());

  return 0;  

}

[[cpp11::register]]
int queuing_off_cpp(
  SEXP model,
  std::string model_class
  ) {

  // Getting the right class
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    ptr->queuing_off();
  
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    ptr->queuing_off();
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    ptr->queuing_off();    

  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    ptr->queuing_off();

  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    ptr->queuing_off();
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());

  return 0;  

}


[[cpp11::register]]
data_frame get_hist_total_cpp(
  SEXP model,
  std::string model_class
) {
  
  // Making some room
  std::vector< int > dates;
  std::vector< std::string > status;
  std::vector< int > counts;
  
  // Getting the right class
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    ptr->get_db().get_hist_total(&dates, &status, &counts);
  
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    ptr->get_db().get_hist_total(&dates, &status, &counts);
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    ptr->get_db().get_hist_total(&dates, &status, &counts);
    
  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    ptr->get_db().get_hist_total(&dates, &status, &counts);
    
  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    ptr->get_db().get_hist_total(&dates, &status, &counts);
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());
  
  
  // Preparing the output
  writable::data_frame res({
    "dates"_nm = dates,
    "status"_nm = status,
    "counts"_nm = counts
  });
  
  
  return res;
}

[[cpp11::register]]
doubles get_transition_probability_cpp(
    SEXP model,
    std::string model_class
) {
  
  // Making some room

  // Getting the right class
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
    
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
    
  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
    
  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    return cpp11::writable::doubles(ptr->get_db().transition_probability(false));
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());
  

}

  [[cpp11::register]]
cpp11::strings get_status_cpp(
    SEXP model,
    std::string model_class
) {
  
  // Making some room
  
  // Getting the right class
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    return cpp11::writable::strings(ptr->get_status());
    
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    return cpp11::writable::strings(ptr->get_status());
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    return cpp11::writable::strings(ptr->get_status());
    
  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    return cpp11::writable::strings(ptr->get_status());
    
  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    return cpp11::writable::strings(ptr->get_status());
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());
  
  
}
  
#define AMap std::unordered_map< std::vector< int >, int, epiworld::vecHasher<int>>
  
  [[cpp11::register]]
cpp11::data_frame get_reproductive_number_cpp(
    SEXP model,
    std::string model_class
) {
  
  // Making some room
  std::vector< int > variant;
  std::vector< int > source;
  std::vector< int > source_exposure_dates;
  std::vector< int > counts;
  
  // Getting the right class
  AMap rn;
  if (model_class == "epiworld_sir")
  {
    
    WrapModel(ModelSIR, ptr)(model);
    rn = ptr->get_db().reproductive_number();
    
  } else if (model_class == "epiworld_sis")
  {
    
    WrapModel(ModelSIS, ptr)(model);
    rn = ptr->get_db().reproductive_number();
    
  } else if (model_class == "epiworld_seir")
  {
    
    WrapModel(ModelSEIR, ptr)(model);
    rn = ptr->get_db().reproductive_number();
    
  } else if (model_class == "epiworld_seirconn")
  {
    
    WrapModel(ModelSEIRCONN, ptr)(model);
    rn = ptr->get_db().reproductive_number();
    
  } else if (model_class == "epiworld_sirconn")
  {
    
    WrapModel(ModelSIRCONN, ptr)(model);
    rn = ptr->get_db().reproductive_number();
    
  } else // Error
    cpp11::stop("Objects of class %s are not supported", model_class.c_str());
  
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

#undef WrapModel
#undef AMap
