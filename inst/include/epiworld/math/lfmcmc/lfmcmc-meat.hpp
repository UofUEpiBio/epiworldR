#ifndef EPIWORLD_LFMCMC_MEAT_HPP
#define EPIWORLD_LFMCMC_MEAT_HPP

#include "lfmcmc-bones.hpp"

/**
 * @brief Proposal function
 * @param params_now Vector where to save the new parameters.
 * @param params_prev Vector of reference parameters.
 * @param m LFMCMC model.
 * @tparam TData 
 */
template<typename TData>
inline void proposal_fun_normal(
    std::vector< epiworld_double >& params_now,
    const std::vector< epiworld_double >& params_prev,
    LFMCMC<TData>* m
) {

    for (size_t p = 0u; p < m->get_n_parameters(); ++p)
        params_now[p] = params_prev[p] + m->rnorm();

    return;
}

/**
 * @brief Factory for a reflective normal kernel
 * 
 * @details Reflective kernel corrects proposals by forcing them to be
 * within prespecified boundaries. 
 * 
 * @tparam TData 
 * @param scale Scale of the normal kernel
 * @param lb Lower bound (applies the same to all parameters)
 * @param ub Upper bound (applies the same to all parameters)
 * @return LFMCMCProposalFun<TData> 
 */
template<typename TData>
inline LFMCMCProposalFun<TData> make_proposal_norm_reflective(
    epiworld_double scale,
    epiworld_double lb,
    epiworld_double ub
) {

    LFMCMCProposalFun<TData> fun =
        [scale,lb,ub](
            std::vector< epiworld_double >& params_now,
            const std::vector< epiworld_double >& params_prev,
            LFMCMC<TData>* m
        ) {

        // Making the proposal
        for (size_t p = 0u; p < m->get_n_parameters(); ++p)
            params_now[p] = params_prev[p] + m->rnorm() * scale;

        // Checking boundaries
        epiworld_double d = ub - lb;
        int odd;
        epiworld_double d_above, d_below;
        for (auto & p : params_now)
        {

            // Correcting if parameter goes above the upper bound
            if (p > ub)
            {
                d_above = p - ub;
                odd     = static_cast<int>(std::floor(d_above / d)) % 2;
                d_above = d_above - std::floor(d_above / d) * d;

                p = (lb + d_above) * odd +
                    (ub - d_above) * (1 - odd);

            // Correcting if parameter goes below upper bound
            } else if (p < lb)
            {
                d_below = lb - p;
                int odd = static_cast<int>(std::floor(d_below / d)) % 2;
                d_below = d_below - std::floor(d_below / d) * d;

                p = (ub - d_below) * odd +
                    (lb + d_below) * (1 - odd);
            }

        }

        #ifdef EPI_DEBUG
        for (auto & p : params_now)
            if (p < lb || p > ub)
                throw std::range_error("The parameter is out of bounds.");
        #endif


        return;

    };

    return fun;
}

/**
 * @brief Uniform proposal kernel
 * 
 * Proposals are made within a radious 1 of the current
 * state of the parameters.
 * 
 * @param params_now Where to write the new parameters
 * @param params_prev Reference parameters
 * @tparam TData 
 * @param m LFMCMC model.
 */
template<typename TData>
inline void proposal_fun_unif(
    std::vector< epiworld_double >& params_now,
    const std::vector< epiworld_double >& params_prev,
    LFMCMC<TData>* m
) {

    for (size_t p = 0u; p < m->get_n_parameters(); ++p)
        params_now[p] = (params_prev[p] + m->runif(-1.0, 1.0));

    return;
}

/**
 * @brief Uses the uniform kernel with euclidean distance
 * 
 * @param stats_now Vector of current statistics based on 
 * simulated data.
 * @param stats_obs Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model.
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_uniform(
    const std::vector< epiworld_double >& stats_now,
    const std::vector< epiworld_double >& stats_obs,
    epiworld_double epsilon,
    LFMCMC<TData>* m
) {

    epiworld_double ans = 0.0;
    for (size_t p = 0u; p < m->get_n_parameters(); ++p)
        ans += std::pow(stats_obs[p] - stats_now[p], 2.0);

    return std::sqrt(ans) < epsilon ? 1.0 : 0.0;

}

#define sqrt2pi() 2.5066282746310002416

/**
 * @brief Gaussian kernel
 * 
 * @tparam TData 
 * @param epsilon 
 * @param m 
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_gaussian(
    const std::vector< epiworld_double >& stats_now,
    const std::vector< epiworld_double >& stats_obs,
    epiworld_double epsilon,
    LFMCMC<TData>* m
) {

    epiworld_double ans = 0.0;
    for (size_t p = 0u; p < m->get_n_parameters(); ++p)
        ans += std::pow(stats_obs[p] - stats_now[p], 2.0);

    return std::exp(
        -.5 * (ans/std::pow(1 + std::pow(epsilon, 2.0)/3.0, 2.0))
        ) / sqrt2pi() ;

}


template<typename TData>
inline void LFMCMC<TData>::set_proposal_fun(LFMCMCProposalFun<TData> fun)
{
    proposal_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_simulation_fun(LFMCMCSimFun<TData> fun)
{
    simulation_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_summary_fun(LFMCMCSummaryFun<TData> fun)
{
    summary_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_kernel_fun(LFMCMCKernelFun<TData> fun)
{
    kernel_fun = fun;
}


template<typename TData>
inline void LFMCMC<TData>::run(
    std::vector< epiworld_double > params_init_,
    size_t n_samples_,
    epiworld_double epsilon_
    )
{

    // Starting timing
    chrono_start();

    // Setting the baseline parameters of the model
    n_samples    = n_samples_;
    epsilon      = epsilon_;
    params_init  = params_init_;
    n_parameters = params_init_.size();

    params_now.resize(n_parameters);
    params_prev.resize(n_parameters);

    if (sampled_data != nullptr)
        sampled_data->resize(n_samples);

    params_prev = params_init;
    params_now  = params_init;

    // Computing the baseline sufficient statistics
    summary_fun(observed_stats, observed_data, this);
    n_statistics = observed_stats.size();

    // Reserving size
    drawn_prob.resize(n_samples);
    sampled_accepted.resize(n_samples, false);
    sampled_stats.resize(n_samples * n_statistics);
    sampled_stats_prob.resize(n_samples);

    accepted_params.resize(n_samples * n_parameters);
    accepted_stats.resize(n_samples * n_statistics);
    accepted_params_prob.resize(n_samples);

    TData data_i = simulation_fun(params_init, this);

    std::vector< epiworld_double > proposed_stats_i;
    summary_fun(proposed_stats_i, data_i, this);
    accepted_params_prob[0u] = kernel_fun(proposed_stats_i, observed_stats, epsilon, this);

    // Recording statistics
    for (size_t i = 0u; i < n_statistics; ++i)
        sampled_stats[i] = proposed_stats_i[i];

    for (size_t k = 0u; k < n_statistics; ++k)
        accepted_params[k] = proposed_stats_i[k];
   
    for (size_t i = 1u; i < n_samples; ++i)
    {
        // Step 1: Generate a proposal and store it in params_now
        proposal_fun(params_now, params_prev, this);

        // Step 2: Using params_now, simulate data
        TData data_i = simulation_fun(params_now, this);

        // Are we storing the data?
        if (sampled_data != nullptr)
            sampled_data->operator[](i) = data_i;

        // Step 3: Generate the summary statistics of the data
        summary_fun(proposed_stats_i, data_i, this);

        // Step 4: Compute the hastings ratio using the kernel function
        epiworld_double hr = kernel_fun(proposed_stats_i, observed_stats, epsilon, this);
        sampled_stats_prob[i] = hr;

        // Storing data
        for (size_t k = 0u; k < n_statistics; ++k)
            sampled_stats[i * n_statistics + k] = proposed_stats_i[k];
        
        // Running Hastings ratio
        epiworld_double r      = runif();
        drawn_prob[i] = r;

        // Step 5: Update if likely
        if (r < std::min(static_cast<epiworld_double>(1.0), hr / accepted_params_prob[i - 1u]))
        {
            accepted_params_prob[i] = hr;
            sampled_accepted[i]     = true;
            
            for (size_t k = 0u; k < n_statistics; ++k)
                accepted_stats[i * n_statistics + k] =
                    proposed_stats_i[k];

            params_prev = params_now;

        } else
        {

            for (size_t k = 0u; k < n_statistics; ++k)
                accepted_stats[i * n_statistics + k] =
                    accepted_stats[(i - 1) * n_statistics + k];

            accepted_params_prob[i] = accepted_params_prob[i - 1u];
        }
            

        for (size_t k = 0u; k < n_parameters; ++k)
            accepted_params[i * n_parameters + k] = params_prev[k];

    }

    // End timing
    chrono_end();

}


template<typename TData>
inline epiworld_double LFMCMC<TData>::runif()
{
    return runifd->operator()(*engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::runif(
    epiworld_double lb,
    epiworld_double ub
)
{
    return runifd->operator()(*engine) * (ub - lb) + lb;
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rnorm()
{
    return rnormd->operator()(*engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rnorm(
    epiworld_double mean,
    epiworld_double sd
    )
{
    return (rnormd->operator()(*engine)) * sd + mean;
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rgamma()
{
    return rgammad->operator()(*engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rgamma(
    epiworld_double alpha,
    epiworld_double beta
    )
{

    auto old_param = rgammad->param();

    rgammad->param(std::gamma_distribution<>::param_type(alpha, beta));

    epiworld_double ans = rgammad->operator()(*engine);

    rgammad->param(old_param);

    return ans;

}

template<typename TData>
inline void LFMCMC<TData>::seed(epiworld_fast_uint s) {

    this->engine->seed(s);

}

template<typename TData>
inline void LFMCMC<TData>::set_rand_engine(std::shared_ptr< std::mt19937 > & eng)
{
    engine = eng;
}

template<typename TData>
inline void LFMCMC<TData>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::make_shared<std::gamma_distribution<>>(alpha,beta);
}

template<typename TData>
inline std::shared_ptr< std::mt19937 > & LFMCMC<TData>::get_rand_endgine()
{
    return engine;
}

// Step 1: Simulate data

// Step 2: Compute the sufficient statistics

// Step 3: Compute the hastings-ratio

// Step 4: Accept/reject, and go back to step 1

#define DURCAST(tunit,txtunit) {\
        elapsed       = std::chrono::duration_cast<std::chrono:: tunit>(\
            time_end - time_start).count(); \
        abbr_unit     = txtunit;}

template<typename TData>
inline void LFMCMC<TData>::get_elapsed(
    std::string unit,
    epiworld_double * last_elapsed,
    std::string * unit_abbr,
    bool print
) {

    // Preparing the result
    epiworld_double elapsed;
    std::string abbr_unit;

    // Figuring out the length
    if (unit == "auto")
    {

        size_t tlength = std::to_string(
            static_cast<int>(floor(time_elapsed.count()))
            ).length();
        
        if (tlength <= 1)
            unit = "nanoseconds";
        else if (tlength <= 3)
            unit = "microseconds";
        else if (tlength <= 6)
            unit = "milliseconds";
        else if (tlength <= 8)
            unit = "seconds";
        else if (tlength <= 9)
            unit = "minutes";
        else 
            unit = "hours";

    }

    if (unit == "nanoseconds")       DURCAST(nanoseconds,"ns")
    else if (unit == "microseconds") DURCAST(microseconds,"\xC2\xB5s")
    else if (unit == "milliseconds") DURCAST(milliseconds,"ms")
    else if (unit == "seconds")      DURCAST(seconds,"s")
    else if (unit == "minutes")      DURCAST(minutes,"m")
    else if (unit == "hours")        DURCAST(hours,"h")
    else
        throw std::range_error("The time unit " + unit + " is not supported.");


    if (last_elapsed != nullptr)
        *last_elapsed = elapsed;
    if (unit_abbr != nullptr)
        *unit_abbr = abbr_unit;

    if (!print)
        return;

    printf_epiworld("Elapsed time : %.2f%s.\n", elapsed, abbr_unit.c_str());
}

#undef DURCAST

#include "lfmcmc-meat-print.hpp"

template<typename TData>
inline void LFMCMC<TData>::chrono_start() {
    time_start = std::chrono::steady_clock::now();
}

template<typename TData>
inline void LFMCMC<TData>::chrono_end() {
    time_end = std::chrono::steady_clock::now();
    time_elapsed += (time_end - time_start);
}

template<typename TData>
inline void LFMCMC<TData>::set_par_names(std::vector< std::string > names)
{

    if (names.size() != n_parameters)
        throw std::length_error("The number of names to add differs from the number of parameters in the model.");

    names_parameters = names;

}
template<typename TData>
inline void LFMCMC<TData>::set_stats_names(std::vector< std::string > names)
{

    if (names.size() != n_statistics)
        throw std::length_error("The number of names to add differs from the number of statistics in the model.");

    names_statistics = names;

}

template<typename TData>
inline std::vector< epiworld_double > LFMCMC<TData>::get_params_mean()
{
    std::vector< epiworld_double > res(this->n_parameters, 0.0);
    
    for (size_t k = 0u; k < n_parameters; ++k)
    {
        for (size_t i = 0u; i < n_samples; ++i)
            res[k] += (this->accepted_params[k + n_parameters * i])/
                static_cast< epiworld_double >(n_samples);
    }

    return res;

}

template<typename TData>
inline std::vector< epiworld_double > LFMCMC<TData>::get_stats_mean()
{
    std::vector< epiworld_double > res(this->n_statistics, 0.0);
    
    for (size_t k = 0u; k < n_statistics; ++k)
    {
        for (size_t i = 0u; i < n_samples; ++i)
            res[k] += (this->accepted_stats[k + n_statistics * i])/
                static_cast< epiworld_double >(n_samples);
    }

    return res;

}

#endif