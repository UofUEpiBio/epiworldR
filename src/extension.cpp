// The build stamp of this epiworldR, which C++ extensions compare with their
// own (see inst/include/epiworldR/extension.h)
#include "cpp11.hpp"
#include "epiworldR/extension.h"

extern "C" const char * epiworldR_build_stamp()
{
    static const std::string stamp = epiworldR::build_stamp();
    return stamp.c_str();
}

[[cpp11::init]]
void register_build_stamp(DllInfo * dll)
{
    R_RegisterCCallable(
        "epiworldR", "build_stamp",
        reinterpret_cast< DL_FUNC >(&epiworldR_build_stamp)
    );
}
