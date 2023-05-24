
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP globalaction_tool_logit_cpp(
    SEXP tool,
    std::vector< int > vars,
    std::vector< double > coefs
) {

    std::vector< size_t > vars_size_t(vars.begin(), vars.end());

    cpp11::external_pointer<GlobalFun<int>> ptr(
        new GlobalFun<int>(globalaction_tool_logit<int>(
            *cpp11::external_pointer<Tool<int>>(tool),
            vars_size_t,
            coefs
        ))
    );
    
    return ptr;

}

[[cpp11::register]]
SEXP globalaction_tool_cpp(
    SEXP tool,
    double prob
) {
    
    cpp11::external_pointer<GlobalFun<int>> ptr(
        new GlobalFun<int>(globalaction_tool<int>(
            *cpp11::external_pointer<Tool<int>>(tool),
            prob
        ))
    );
    
    return ptr;

}

[[cpp11::register]]
SEXP add_global_action_cpp(
    SEXP model,
    SEXP action,
    int date
) {
  
  external_pointer<Model<int>> modelptr(model);
  external_pointer<GlobalFun<int>> actionptr(action);
  
  modelptr->add_global_action(*actionptr, date);
  
  return model;
  
}
