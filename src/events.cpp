
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

    GlobalFun event(epimodels::globalevent_tool_logit(
        *cpp11::external_pointer<Tool>(tool),
        vars_size_t,
        coefs
    ));

    cpp11::external_pointer<GlobalEvent> ptr(
        new GlobalEvent(event, name, day)
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

    GlobalFun event(epimodels::globalevent_tool(
        *cpp11::external_pointer<Tool>(tool),
        prob
    ));

    cpp11::external_pointer<GlobalEvent> ptr(
        new GlobalEvent(event, name, day)
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

    GlobalFun event(epimodels::globalevent_set_param(
        param,
        value
    ));

    cpp11::external_pointer<GlobalEvent> ptr(
        new GlobalEvent(event, name, day)
    );

    return ptr;

}

[[cpp11::register]]
SEXP print_globalevent_cpp(
    SEXP event
) {

  external_pointer<GlobalEvent> eventptr(event);

  eventptr->print();

  return event;

}


[[cpp11::register]]
SEXP add_globalevent_cpp(
    SEXP model,
    SEXP event
) {

  external_pointer<Model> modelptr(model);
  external_pointer<GlobalEvent> eventptr(event);

  modelptr->add_globalevent(*eventptr);

  return model;

}

[[cpp11::register]]
SEXP rm_globalevent_cpp(
    SEXP model,
    std::string name
) {

  external_pointer<Model> modelptr(model);

  modelptr->rm_globalevent(name);

  return model;

}

[[cpp11::register]]
SEXP globalevent_fun_cpp(
    cpp11::function fun,
    std::string name,
    int day
    ) {

  // NOTE: This lambda calls the R API. R is NOT thread-safe and must only be
  // called from the main thread. When used inside run_multiple() with
  // nthreads > 1, this lambda is invoked from OpenMP worker threads, which
  // is undefined behavior. Users must set nthreads = 1 when using
  // globalevent_fun() with run_multiple().
  GlobalFun fun_call = [fun](Model * model) -> void {

    cpp11::external_pointer<Model> modelptr(model, false);

    sexp modelptrs(modelptr);
    modelptrs.attr("class") = "epiworld_model";

    fun(modelptr);

    return;

  };

  return external_pointer<GlobalEvent>(
    new GlobalEvent(fun_call, name, day)
  );


}
