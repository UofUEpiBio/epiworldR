#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/r_vector.hpp"
#include "cpp11/sexp.hpp"
#include <iostream>

#include "epiworld-common.h"

using namespace epiworld;

#define TData_default std::vector< int >

#define WrapLFMCMC(a) \
    cpp11::external_pointer<LFMCMC<TData_default>> (a)

// LFMCMC definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/math/lfmcmc

[[cpp11::register]]
SEXP LFMCMC_cpp(
    SEXP model
) {
    WrapLFMCMC(lfmcmc_ptr)(
        new LFMCMC<TData_default>()
    );

    cpp11::external_pointer<Model<int>> modelptr(model);
    lfmcmc_ptr->set_rand_engine(modelptr->get_rand_endgine());

    return lfmcmc_ptr;
}

[[cpp11::register]]
SEXP run_lfmcmc_cpp(
    SEXP lfmcmc,
    std::vector<epiworld_double> params_init_,
    size_t n_samples_,
    epiworld_double epsilon_,
    int seed
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->run(params_init_, n_samples_, epsilon_, seed);
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
    cpp11::function fun
) {

    cpp11::stop("Un implemented");

    return lfmcmc;
}

// Use proposal function defined in epiworld
[[cpp11::register]]
SEXP use_proposal_norm_reflective_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_proposal_fun(make_proposal_norm_reflective<std::vector<int>>(.5, 0, 1));
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_simulation_fun_cpp(
    SEXP lfmcmc,
    cpp11::function fun
) {

    LFMCMCSimFun<TData_default> fun_call = [fun](
        const std::vector<epiworld_double>& params,
        LFMCMC<TData_default>*
        ) -> TData_default {

        auto params_doubles = cpp11::doubles(params);

        return cpp11::as_cpp<TData_default>(
            cpp11::integers(fun(params_doubles))
            );
    };

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);

    lfmcmc_ptr->set_simulation_fun(fun_call);

    return lfmcmc;
}

[[cpp11::register]]
SEXP set_summary_fun_cpp(
    SEXP lfmcmc,
    cpp11::function fun
) {

    LFMCMCSummaryFun<TData_default> fun_call = [fun](
        std::vector< epiworld_double >& res,
        const TData_default& dat,
        LFMCMC<TData_default>*
        ) -> void {

        if (res.size() == 0u)
            res.resize(dat.size());

        auto dat_int = cpp11::integers(dat);
        auto res_tmp = cpp11::integers(fun(dat_int));

        std::copy(res_tmp.begin(), res_tmp.end(), res.begin());

        return;
    };

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);

    lfmcmc_ptr->set_summary_fun(fun_call);

    return lfmcmc;
}

[[cpp11::register]]
SEXP set_kernel_fun_cpp(
    SEXP lfmcmc,
    cpp11::function fun
) {
    cpp11::stop("Unimplemented");

    return lfmcmc;
}

// Use kernel function defined in epiworld
[[cpp11::register]]
SEXP use_kernel_fun_gaussian_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_kernel_fun(kernel_fun_gaussian<std::vector<int>>);
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
