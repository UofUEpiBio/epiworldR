#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/r_vector.hpp"

#include "epiworld-common.h"

using namespace epiworld;

#define TData_default std::vector< int >
// #define LFMCMCSimFun_default std::function<TData_default(const std::vector< epiworld_double >&,epiworld::LFMCMC<TData_default>*)>
// #define LFMCMCSummaryFun_default std::function<void(std::vector< epiworld_double >&,const TData_default&,epiworld::LFMCMC<TData_default>*)>
// #define LFMCMCProposalFun_default std::function<void(std::vector< epiworld_double >&,const std::vector< epiworld_double >&,epiworld::LFMCMC<TData_default>*)>
// #define LFMCMCKernelFun_default std::function<epiworld_double(const std::vector< epiworld_double >&,const std::vector< epiworld_double >&,epiworld_double,epiworld::LFMCMC<TData_default>*)>

#define WrapLFMCMC(a) \
  cpp11::external_pointer<LFMCMC<TData_default>> (a)

// LFMCMC definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/math/lfmcmc

[[cpp11::register]]
SEXP LFMCMC_cpp() {
    WrapLFMCMC(lfmcmc_ptr)(
        new LFMCMC<TData_default>()
    );

    return lfmcmc_ptr;
}

[[cpp11::register]]
SEXP run_lfmcmc_cpp(
    SEXP lfmcmc,
    std::vector<epiworld_double> params_init_,
    size_t n_samples_,
    epiworld_double epsilon_
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->run(params_init_, n_samples_, epsilon_);
    return lfmcmc;
}

// observed_data_ should be of type TData
[[cpp11::register]]
SEXP set_observed_data_cpp(
    SEXP lfmcmc,
    std::vector< int > observed_data_
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_observed_data(observed_data_);
    return lfmcmc;
}

// [[cpp11::register]]
// SEXP set_proposal_fun_cpp(
//     SEXP lfmcmc,
//     LFMCMCProposalFun_default<TData_default> fun
// ) {
//     WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
//     lfmcmc_ptr->set_proposal_fun(fun);
//     return lfmcmc;
// }

// [[cpp11::register]]
// SEXP set_simulation_fun_cpp(
//     SEXP lfmcmc,
//     LFMCMCSimFun_default<TData_default> fun
// ) {
//     WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
//     lfmcmc_ptr->set_simulation_fun(fun);
//     return lfmcmc;
// }

// [[cpp11::register]]
// SEXP set_summary_fun_cpp(
//     SEXP lfmcmc,
//     LFMCMCSummaryFun_default<TData_default> fun
// ) {
//     WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
//     lfmcmc_ptr->set_summary_fun(fun);
//     return lfmcmc;
// }

// [[cpp11::register]]
// SEXP set_kernel_fun_cpp(
//     SEXP lfmcmc,
//     LFMCMCKernelFun_default<TData_default> fun
// ) {
//     WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
//     lfmcmc_ptr->set_kernel_fun(fun);
//     return lfmcmc;
// }

#undef WrapLFMCMC
