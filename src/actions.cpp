
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
    std::vector< double > coefs,
    std::string name,
    int day
) {

    std::vector< size_t > vars_size_t(vars.begin(), vars.end());

    GlobalFun<int> action(epimodels::globalaction_tool_logit<int>(
        *cpp11::external_pointer<Tool<int>>(tool),
        vars_size_t,
        coefs
    ));

    cpp11::external_pointer<GlobalAction<int>> ptr(
        new GlobalAction<int>(action, name, day)
    );
    
    return ptr;

}

[[cpp11::register]]
SEXP globalaction_tool_cpp(
    SEXP tool,
    double prob,
    std::string name,
    int day
) {

    GlobalFun<int> action(epimodels::globalaction_tool<int>(
        *cpp11::external_pointer<Tool<int>>(tool),
        prob
    ));
    
    cpp11::external_pointer<GlobalAction<int>> ptr(
        new GlobalAction<int>(action, name, day)
    );
    
    return ptr;

}

[[cpp11::register]]
SEXP globalaction_set_param_cpp(
    std::string param,
    double value,
    std::string name,
    int day
) {

    GlobalFun<int> action(epimodels::globalaction_set_param<int>(
        param,
        value
    ));
    
    cpp11::external_pointer<GlobalAction<int>> ptr(
        new GlobalAction<int>(action, name, day)
    );
    
    return ptr;

}

[[cpp11::register]]
SEXP print_global_action_cpp(
    SEXP action
) {
  
  external_pointer<GlobalAction<int>> actionptr(action);
  
  actionptr->print();
  
  return action;
  
}


[[cpp11::register]]
SEXP add_global_action_cpp(
    SEXP model,
    SEXP action
) {
  
  external_pointer<Model<int>> modelptr(model);
  external_pointer<GlobalAction<int>> actionptr(action);
  
  modelptr->add_global_action(*actionptr);
  
  return model;
  
}

[[cpp11::register]]
SEXP rm_global_action_cpp(
    SEXP model,
    std::string name
) {
  
  external_pointer<Model<int>> modelptr(model);
  
  modelptr->rm_global_action(name);
  
  return model;
  
}

