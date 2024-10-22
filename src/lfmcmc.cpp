#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/r_vector.hpp"

#include "epiworld-common.h"

using namespace epiworld;

#define TData_default std::vector< int >

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

[[cpp11::register]]
SEXP set_proposal_fun_cpp(
    SEXP lfmcmc,
    SEXP fun
) {
    cpp11::external_pointer<LFMCMCProposalFun<TData_default>> fun_ptr(fun);
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_proposal_fun(*fun_ptr);
    return lfmcmc;
}

[[cpp11::register]]
SEXP create_LFMCMCSimFun_cpp(
    cpp11::function fun
    ) {

    LFMCMCSimFun<TData_default> fun_call = [fun](std::vector<epiworld_double>& params, LFMCMC<TData_default>* model) {
        WrapLFMCMC(lfmcmc_ptr)(model);
        SEXP res = fun(params, lfmcmc_ptr);
        cpp11::external_pointer<TData_default> res_vec(res);
        return *res_vec;
    };

    return cpp11::external_pointer<LFMCMCSimFun<TData_default>>(
        new LFMCMCSimFun<TData_default>(fun_call)
    );
}

[[cpp11::register]]
SEXP set_simulation_fun_cpp(
    SEXP lfmcmc,
    SEXP fun
) {
    cpp11::external_pointer<LFMCMCSimFun<TData_default>> fun_ptr(fun);
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_simulation_fun(*fun_ptr);
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_summary_fun_cpp(
    SEXP lfmcmc,
    SEXP fun
) {
    cpp11::external_pointer<LFMCMCSummaryFun<TData_default>> fun_ptr(fun);
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_summary_fun(*fun_ptr);
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_kernel_fun_cpp(
    SEXP lfmcmc,
    SEXP fun
) {
    cpp11::external_pointer<LFMCMCKernelFun<TData_default>> fun_ptr(fun);
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_kernel_fun(*fun_ptr);
    return lfmcmc;
}

// s should be of type epiworld_fast_uint
[[cpp11::register]]
SEXP seed_lfmcmc_cpp(
    SEXP lfmcmc,
    unsigned long long int s
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->seed(s);
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_par_names_cpp(
    SEXP lfmcmc,
    std::vector< std::string > names
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_par_names(names);
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_stats_names_cpp(
    SEXP lfmcmc,
    std::vector< std::string > names
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_stats_names(names);
    return lfmcmc;
}

[[cpp11::register]]
SEXP print_lfmcmc_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->print();
    return lfmcmc;
}

#undef WrapLFMCMC
