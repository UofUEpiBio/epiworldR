#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"

#include "epiworld-common.h"

// LFMCMC definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/math/lfmcmc

[[cpp11::register]]
SEXP LFMCMC_cpp() {
    cpp11::external_pointer<epiworld::LFMCMC<>> lfmcmc_ptr(
        new epiworld::LFMCMC<>()
    );

    return lfmcmc_ptr;
}

