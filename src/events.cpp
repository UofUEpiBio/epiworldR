
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/data_frame.hpp"
#include "epiworld-common.h"
#include "cpp11/sexp.hpp"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP globalevent_tool_logit_cpp(
    SEXP tool,
    std::vector< int > vars,
    std::vector< double > coefs,
    std::string name,
    int day
) {

    std::vector< size_t > vars_size_t(vars.begin(), vars.end());

    GlobalFun<int> event(epimodels::globalevent_tool_logit<int>(
        *cpp11::external_pointer<Tool<int>>(tool),
        vars_size_t,
        coefs
    ));

    cpp11::external_pointer<GlobalEvent<int>> ptr(
        new GlobalEvent<int>(event, name, day)
    );

    return ptr;

}

[[cpp11::register]]
SEXP globalevent_tool_cpp(
    SEXP tool,
    double prob,
    std::string name,
    int day
) {

    GlobalFun<int> event(epimodels::globalevent_tool<int>(
        *cpp11::external_pointer<Tool<int>>(tool),
        prob
    ));

    cpp11::external_pointer<GlobalEvent<int>> ptr(
        new GlobalEvent<int>(event, name, day)
    );

    return ptr;

}

[[cpp11::register]]
SEXP globalevent_set_param_cpp(
    std::string param,
    double value,
    std::string name,
    int day
) {

    GlobalFun<int> event(epimodels::globalevent_set_param<int>(
        param,
        value
    ));

    cpp11::external_pointer<GlobalEvent<int>> ptr(
        new GlobalEvent<int>(event, name, day)
    );

    return ptr;

}

[[cpp11::register]]
SEXP print_globalevent_cpp(
    SEXP event
) {

  external_pointer<GlobalEvent<int>> eventptr(event);

  eventptr->print();

  return event;

}


[[cpp11::register]]
SEXP add_globalevent_cpp(
    SEXP model,
    SEXP event
) {

  external_pointer<Model<int>> modelptr(model);
  external_pointer<GlobalEvent<int>> eventptr(event);

  modelptr->add_globalevent(*eventptr);

  return model;

}

[[cpp11::register]]
SEXP rm_globalevent_cpp(
    SEXP model,
    std::string name
) {

  external_pointer<Model<int>> modelptr(model);

  modelptr->rm_globalevent(name);

  return model;

}

[[cpp11::register]]
SEXP globalevent_fun_cpp(
    cpp11::function fun,
    std::string name,
    int day
    ) {

  GlobalFun<int> fun_call = [fun](Model<int> * model) -> void {

    cpp11::external_pointer<Model<int>> modelptr(model, false);

    sexp modelptrs(modelptr);
    modelptrs.attr("class") = "epiworld_model";

    fun(modelptr);

    return;

  };

  return external_pointer<GlobalEvent<int>>(
    new GlobalEvent<int>(fun_call, name, day)
  );


}
