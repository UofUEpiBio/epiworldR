#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"

using namespace epiworld;

#define WrapModelDiagram(a) \
  cpp11::external_pointer<ModelDiagram> (a)

// ModelDiagram definitions:
// https://github.com/UofUEpiBio/epiworld/blob/master/include/epiworld/modeldiagram-bones.hpp

[[cpp11::register]]
SEXP ModelDiagram_cpp() {
    WrapModelDiagram(ptr)(
        new ModelDiagram()
    );

    return ptr;
}

[[cpp11::register]]
SEXP draw_from_data_cpp(
    SEXP model_diagram,
    const std::vector< std::string > & states,
    const std::vector< epiworld_double > & tprob,
    const std::string & fn_output,
    bool self
) {
    WrapModelDiagram(diagram_ptr)(model_diagram);
    diagram_ptr->draw_from_data(
        states,
        tprob,
        fn_output,
        self
    );
    return model_diagram;
}

[[cpp11::register]]
SEXP draw_from_file_cpp(
    SEXP model_diagram,
    const std::string & fn_transition,
    const std::string & fn_output,
    bool self
) {
    WrapModelDiagram(diagram_ptr)(model_diagram);
    diagram_ptr->draw_from_file(
        fn_transition,
        fn_output,
        self
    );
    return model_diagram;
}

[[cpp11::register]]
SEXP draw_from_files_cpp(
    SEXP model_diagram,
    const std::vector< std::string > & fns_transition,
    const std::string & fn_output,
    bool self
) {
    WrapModelDiagram(diagram_ptr)(model_diagram);
    diagram_ptr->draw_from_files(
        fns_transition,
        fn_output,
        self
    );
    return model_diagram;
}
