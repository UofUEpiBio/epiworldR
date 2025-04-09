#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"

using namespace epiworld;

// ModelDiagram definitions:
// https://github.com/UofUEpiBio/epiworld/blob/master/include/epiworld/modeldiagram-bones.hpp

[[cpp11::register]]
void draw_from_data_cpp(
    const std::vector< std::string > & states,
    const std::vector< epiworld_double > & tprob,
    const std::string & fn_output,
    bool self
) {
    epiworld::ModelDiagram diagram;
    diagram.draw_from_data(
        states,
        tprob,
        fn_output,
        self
    );
}

[[cpp11::register]]
void draw_from_file_cpp(
    const std::string & fn_transition,
    const std::string & fn_output,
    bool self
) {
    epiworld::ModelDiagram diagram;
    diagram.draw_from_file(
        fn_transition,
        fn_output,
        self
    );
}

[[cpp11::register]]
void draw_from_files_cpp(
    const std::vector< std::string > & fns_transition,
    const std::string & fn_output,
    bool self
) {
    epiworld::ModelDiagram diagram;
    diagram.draw_from_files(
        fns_transition,
        fn_output,
        self
    );
}
