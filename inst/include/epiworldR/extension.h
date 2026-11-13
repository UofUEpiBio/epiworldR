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
 * The configuration below is epiworldR's (`src/Makevars`): the epiworld
 * classes must be laid out the same in the extension and in epiworldR.
 * `as_model()` checks it, and stops with a message if the extension was built
 * against other headers or with another configuration.
 */

#include <string>
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
#include "epiworld-commit.h"

namespace epiworldR {

/**
 * @brief What an extension must share with epiworldR: the epiworld commit
 * of the headers, the size of `epiworld_double`, the debug switches that add
 * members, and the sizes of the classes an extension uses.
 * @details `static` (one copy per translation unit): with external linkage,
 * the dynamic linker could resolve an extension's call to epiworldR's own
 * copy, which reports epiworldR's stamp, and a mismatch would go unseen.
 */
static inline std::string build_stamp()
{

    auto size = [](const char * name, size_t n) {
        return std::string("; ") + name + " " + std::to_string(n);
    };

    return std::string("epiworld ") + EPIWORLDR_EPIWORLD_COMMIT +
        size("epiworld_double", sizeof(epiworld_double)) +
        size("Model", sizeof(epiworld::Model)) +
        size("Agent", sizeof(epiworld::Agent)) +
        size("Virus", sizeof(epiworld::Virus)) +
        size("Tool", sizeof(epiworld::Tool)) +
        size("GlobalEvent", sizeof(epiworld::GlobalEvent)) +
        size("DataBase", sizeof(epiworld::DataBase))
        #ifdef EPI_DEBUG
        + "; EPI_DEBUG"
        #endif
        #ifdef EPI_DEBUG_VIRUS
        + "; EPI_DEBUG_VIRUS"
        #endif
        ;

}

/**
 * @brief The model behind an epiworldR model (`Model()`, `make_model()`, a
 * named model, ...).
 * @details Everything the C++ API can do to a model is then available: add
 * global events, tools, viruses with their natural histories, layers,
 * parameters... What is added belongs to the model, as if added from R.
 * Stops with an error if `model` is not an epiworldR model, or if this code
 * was built against other epiworld headers or another configuration than
 * the installed epiworldR (rebuild it). `static`, as `build_stamp()`.
 */
static inline epiworld::Model & as_model(SEXP model)
{

    using stamp_fun = const char * (*)();
    static const stamp_fun installed = reinterpret_cast< stamp_fun >(
        R_GetCCallable("epiworldR", "build_stamp")
    );
    static const std::string here = build_stamp();

    const std::string theirs = installed();
    if (here != theirs)
        cpp11::stop(
            "This C++ code was built for another epiworldR. Rebuild it against "
            "the installed epiworldR.\n  This code: %s\n  epiworldR: %s",
            here.c_str(), theirs.c_str()
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
