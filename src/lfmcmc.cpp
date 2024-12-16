#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "cpp11/r_vector.hpp"
#include "cpp11/sexp.hpp"
#include "cpp11/doubles.hpp"
#include <iostream>

#include "epiworld-common.h"

using namespace epiworld;

#define TData_default std::vector< double >

#define WrapLFMCMC(a) \
    cpp11::external_pointer<LFMCMC<TData_default>> (a)

inline cpp11::sexp lfmcmc_as_sexp(LFMCMC<TData_default> * lfmcmc) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc, false);
    cpp11::sexp lfmcmc_ptr_s(lfmcmc_ptr);
    lfmcmc_ptr_s.attr("class") = "epiworld_lfmcmc";
    return lfmcmc_ptr_s;
}

// LFMCMC definitions:
// https://github.com/UofUEpiBio/epiworld/tree/master/include/epiworld/math/lfmcmc

// *************************************
// LFMCMC Function
// *************************************

[[cpp11::register]]
SEXP LFMCMC_cpp(
    SEXP model
) {
    WrapLFMCMC(lfmcmc_ptr)(
        new LFMCMC<TData_default>()
    );

    if (Rf_inherits(model, "epiworld_model")) {
        cpp11::external_pointer<Model<int>> modelptr(model);
        lfmcmc_ptr->set_rand_engine(modelptr->get_rand_endgine());
    } else {
        auto new_ptr = std::make_shared<std::mt19937>(std::mt19937());
        lfmcmc_ptr->set_rand_engine(new_ptr);
    }

    return lfmcmc_ptr;
}

// *************************************
// LFMCMC Run Function
// *************************************

[[cpp11::register]]
SEXP run_lfmcmc_cpp(
    SEXP lfmcmc,
    std::vector<double> params_init_,
    size_t n_samples_,
    double epsilon_,
    int seed
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->run(params_init_, n_samples_, epsilon_, seed);
    return lfmcmc;
}

// *************************************
// LFMCMC Setup Functions
// *************************************

// observed_data_ should be of type TData
[[cpp11::register]]
SEXP set_observed_data_cpp(
    SEXP lfmcmc,
    std::vector< double > observed_data_
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

    LFMCMCProposalFun<TData_default> fun_call = [fun](
        std::vector< epiworld_double >& new_params,
        const std::vector< epiworld_double >& old_params,
        LFMCMC<TData_default> * lfmcmc_obj
        ) -> void {

        auto params_doubles = cpp11::doubles(old_params);
        auto lfmcmc_ptr = lfmcmc_as_sexp(lfmcmc_obj);
        auto res_tmp = cpp11::doubles(fun(params_doubles, lfmcmc_ptr));

        std::copy(
            res_tmp.begin(),
            res_tmp.end(),
            new_params.begin()
            );

        return;
    };

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);

    lfmcmc_ptr->set_proposal_fun(fun_call);

    return lfmcmc;
}

// Use proposal function defined in epiworld
[[cpp11::register]]
SEXP use_proposal_norm_reflective_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_proposal_fun(make_proposal_norm_reflective<std::vector<double>>(.5, 0, 1));
    return lfmcmc;
}

[[cpp11::register]]
SEXP set_simulation_fun_cpp(
    SEXP lfmcmc,
    cpp11::function fun
) {

    LFMCMCSimFun<TData_default> fun_call = [fun](
        const std::vector<epiworld_double>& params,
        LFMCMC<TData_default> * lfmcmc_obj
        ) -> TData_default {

        auto params_doubles = cpp11::doubles(params);
        auto lfmcmc_ptr = lfmcmc_as_sexp(lfmcmc_obj);
        return cpp11::as_cpp<TData_default>(
            cpp11::doubles(fun(params_doubles, lfmcmc_ptr))
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
        LFMCMC<TData_default>* lfmcmc_obj
        ) -> void {

        auto dat_int = cpp11::doubles(dat);
        auto lfmcmc_ptr = lfmcmc_as_sexp(lfmcmc_obj);
        auto res_tmp = cpp11::doubles(fun(dat_int, lfmcmc_ptr));


        if (res.size() == 0u)
            res.resize(res_tmp.size());

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

    LFMCMCKernelFun<TData_default> fun_call = [fun](
        const std::vector< epiworld_double >& simulated_stats,
        const std::vector< epiworld_double >& observed_stats,
        epiworld_double epsilon,
        LFMCMC<TData_default> * lfmcmc_obj
        ) -> epiworld_double {

        auto sim_stats_doubles = cpp11::doubles(simulated_stats);
        auto obs_stats_doubles = cpp11::doubles(observed_stats);

        auto lfmcmc_ptr = lfmcmc_as_sexp(lfmcmc_obj);

        return cpp11::as_cpp<epiworld_double>(
            fun(sim_stats_doubles, obs_stats_doubles, epsilon, lfmcmc_ptr)
            );
    };

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);

    lfmcmc_ptr->set_kernel_fun(fun_call);

    return lfmcmc;
}

// Use kernel function defined in epiworld
[[cpp11::register]]
SEXP use_kernel_fun_gaussian_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_kernel_fun(kernel_fun_gaussian<std::vector<double>>);
    return lfmcmc;
}

// *************************************
// LFMCMC Getters - Params and Stats
// *************************************

[[cpp11::register]]
cpp11::writable::doubles get_mean_params_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_mean_params());
}

[[cpp11::register]]
cpp11::writable::doubles get_mean_stats_cpp(
    SEXP lfmcmc
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_mean_stats());
}

[[cpp11::register]]
cpp11::writable::doubles get_initial_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_initial_params());

}

[[cpp11::register]]
cpp11::writable::doubles get_current_proposed_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_current_proposed_params());

}

[[cpp11::register]]
cpp11::writable::doubles get_current_accepted_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_current_accepted_params());

}

[[cpp11::register]]
cpp11::writable::doubles get_current_proposed_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_current_proposed_stats());

}

[[cpp11::register]]
cpp11::writable::doubles get_current_accepted_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_current_accepted_stats());

}

[[cpp11::register]]
cpp11::writable::doubles get_observed_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_observed_stats());

}

[[cpp11::register]]
SEXP get_all_sample_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_sample_params());

}

[[cpp11::register]]
SEXP get_all_sample_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_sample_stats());

}

[[cpp11::register]]
cpp11::writable::doubles get_all_sample_acceptance_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_sample_acceptance());

}

[[cpp11::register]]
cpp11::writable::doubles get_all_sample_drawn_prob_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_sample_drawn_prob());

}

[[cpp11::register]]
cpp11::writable::doubles get_all_sample_kernel_scores_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_sample_kernel_scores());

}

[[cpp11::register]]
SEXP get_all_accepted_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_accepted_params());

}

[[cpp11::register]]
SEXP get_all_accepted_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_accepted_stats());

}

[[cpp11::register]]
SEXP get_all_accepted_kernel_scores_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return cpp11::doubles(lfmcmc_ptr->get_all_accepted_kernel_scores());

}

// *************************************
// LFMCMC Getters - N values
// *************************************

[[cpp11::register]]
int get_n_samples_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return static_cast<int>(lfmcmc_ptr->get_n_samples());

}

[[cpp11::register]]
int get_n_stats_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return static_cast<int>(lfmcmc_ptr->get_n_stats());

}

[[cpp11::register]]
int get_n_params_cpp(SEXP lfmcmc) {

    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    return static_cast<int>(lfmcmc_ptr->get_n_params());

}

// *************************************
// LFMCMC Printing functions
// *************************************

[[cpp11::register]]
SEXP verbose_off_lfmcmc_cpp(SEXP lfmcmc) {

  WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
  lfmcmc_ptr->verbose_off();
  return lfmcmc;

}

[[cpp11::register]]
SEXP verbose_on_lfmcmc_cpp(SEXP lfmcmc) {

  WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
  lfmcmc_ptr->verbose_on();
  return lfmcmc;

}

[[cpp11::register]]
SEXP set_params_names_cpp(
    SEXP lfmcmc,
    std::vector< std::string > names
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->set_params_names(names);
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
    SEXP lfmcmc,
    int burnin
) {
    WrapLFMCMC(lfmcmc_ptr)(lfmcmc);
    lfmcmc_ptr->print(static_cast<size_t>(burnin));
    return lfmcmc;
}

#undef WrapLFMCMC
