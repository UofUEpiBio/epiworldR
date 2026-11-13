// The major version of epiworld in this epiworldR, which C++ extensions
// compare with their own (see inst/include/epiworldR/extension.h)
#include "cpp11.hpp"
#include "epiworldR/extension.h"

extern "C" int epiworldR_epiworld_version_major()
{
    return EPIWORLD_VERSION_MAJOR;
}

[[cpp11::init]]
void register_epiworld_version_major(DllInfo * dll)
{
    R_RegisterCCallable(
        "epiworldR", "epiworld_version_major",
        reinterpret_cast< DL_FUNC >(&epiworldR_epiworld_version_major)
    );
}
