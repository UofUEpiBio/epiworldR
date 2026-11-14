#ifndef EPIWORLDR_EXTENSION_H
#define EPIWORLDR_EXTENSION_H

/**
 * @file extension.h
 * @brief Extending epiworldR in C++: reach a model made in R and use the
 * whole epiworld C++ API on it.
 *
 * @details Include this header (and nothing else of epiworld) in C++ that
 * works on epiworldR's models, compiled in a session with
 * `epiworldR::cpp_source_epiworld()` or in a package with
 * `LinkingTo: cpp11, epiworldR`:
 *
 * ```cpp
 * #include <epiworldR/extension.h>
 *
 * [[cpp11::register]]
 * void add_my_event(SEXP model) {
 *     epiworld::Model & m = epiworldR::as_model(model);
 *     m.add_globalevent(...);
 * }
 * ```
 *
 * The configuration below is epiworldR's (`src/Makevars`); do not change
 * it. `as_model()` checks that the extension and the installed epiworldR use
 * the same major version of epiworld.
 */

#include <cpp11.hpp>
#include <R_ext/Rdynload.h>

#ifndef printf_epiworld
    #define printf_epiworld Rprintf
#endif
#ifndef epiworld_double
    #define epiworld_double double
#endif
#ifndef EPI_CHECK_USER_INTERRUPT
    #define EPI_CHECK_USER_INTERRUPT(a) cpp11::check_user_interrupt()
#endif

#include "../epiworld/epiworld.hpp"

namespace epiworldR {

/**
 * @brief The model behind an epiworldR model (`Model()`, `make_model()`, a
 * named model, ...).
 * @details Everything the C++ API can do to a model is then available: add
 * global events, tools, viruses with their natural histories, layers,
 * parameters... What is added belongs to the model, as if added from R.
 * Stops with an error if `model` is not an epiworldR model, or if this code
 * was built against another major version of epiworld than the installed
 * epiworldR's (rebuild it).
 * `static` (one copy per translation unit): with external linkage, the
 * dynamic linker could resolve an extension's call to epiworldR's own copy,
 * which compares epiworldR's version with itself.
 */
static inline epiworld::Model & as_model(SEXP model)
{

    using version_fun = int (*)();
    static const int installed = reinterpret_cast< version_fun >(
        R_GetCCallable("epiworldR", "epiworld_version_major")
    )();

    if (installed != EPIWORLD_VERSION_MAJOR)
        cpp11::stop(
            "This C++ code was built for epiworld %d, but epiworldR uses "
            "epiworld %d. Rebuild it against the installed epiworldR.",
            EPIWORLD_VERSION_MAJOR, installed
        );

    if (!Rf_inherits(model, "epiworld_model"))
        cpp11::stop("`model` is not an epiworldR model.");

    cpp11::external_pointer< epiworld::Model > ptr(model);
    if (ptr.get() == nullptr)
        cpp11::stop(
            "The model is empty (an external pointer does not survive "
            "saving and loading the R object): build it again."
        );

    return *ptr;

}

}

#endif
