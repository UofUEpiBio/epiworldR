#ifndef EPIWORLD_LFMCMC_MEAT_HPP
#define EPIWORLD_LFMCMC_MEAT_HPP

#include "lfmcmc-bones.hpp"

/**
 * @brief Proposal function
 * @param new_params Vector where to save the new parameters.
 * @param old_params Vector of reference parameters.
 * @param m LFMCMC model.
 * @tparam TData 
 */
template<typename TData>
inline void proposal_fun_normal(
    std::vector< epiworld_double >& new_params,
    const std::vector< epiworld_double >& old_params,
    LFMCMC<TData>* m
) {

    for (size_t p = 0u; p < m->get_n_params(); ++p)
        new_params[p] = old_params[p] + m->rnorm();

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
            std::vector< epiworld_double >& new_params,
            const std::vector< epiworld_double >& old_params,
            LFMCMC<TData>* m
        ) {

        // Making the proposal
        for (size_t p = 0u; p < m->get_n_params(); ++p)
            new_params[p] = old_params[p] + m->rnorm() * scale;

        // Checking boundaries
        epiworld_double d = ub - lb;
        int odd;
        epiworld_double d_above, d_below;
        for (auto & p : new_params)
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
        for (auto & p : new_params)
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
 * @param new_params Where to write the new parameters
 * @param old_params Reference parameters
 * @tparam TData 
 * @param m LFMCMC model.
 */
template<typename TData>
inline void proposal_fun_unif(
    std::vector< epiworld_double >& new_params,
    const std::vector< epiworld_double >& old_params,
    LFMCMC<TData>* m
) {

    for (size_t p = 0u; p < m->get_n_params(); ++p)
        new_params[p] = (old_params[p] + m->runif(-1.0, 1.0));

    return;
}

/**
 * @brief Uses the uniform kernel with euclidean distance
 * 
 * @param simulated_stats Vector of statistics based on 
 * simulated data
 * @param observed_stats Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_uniform(
    const std::vector< epiworld_double >& simulated_stats,
    const std::vector< epiworld_double >& observed_stats,
    epiworld_double epsilon,
    LFMCMC<TData>* m
) {

    epiworld_double ans = 0.0;
    for (size_t p = 0u; p < m->get_n_params(); ++p)
        ans += std::pow(observed_stats[p] - simulated_stats[p], 2.0);

    return std::sqrt(ans) < epsilon ? 1.0 : 0.0;

}

#define sqrt2pi() 2.5066282746310002416

/**
 * @brief Gaussian kernel
 * 
 * @tparam TData 
 * @param simulated_stats Vector of statistics based on 
 * simulated data
 * @param observed_stats Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_gaussian(
    const std::vector< epiworld_double >& simulated_stats,
    const std::vector< epiworld_double >& observed_stats,
    epiworld_double epsilon,
    LFMCMC<TData>* m
) {

    epiworld_double ans = 0.0;
    for (size_t p = 0u; p < m->get_n_params(); ++p)
        ans += std::pow(observed_stats[p] - simulated_stats[p], 2.0);

    return std::exp(
        -.5 * (ans/std::pow(1 + std::pow(epsilon, 2.0)/3.0, 2.0))
        ) / sqrt2pi() ;

}


template<typename TData>
inline void LFMCMC<TData>::set_proposal_fun(LFMCMCProposalFun<TData> fun)
{
    m_proposal_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_simulation_fun(LFMCMCSimFun<TData> fun)
{
    m_simulation_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_summary_fun(LFMCMCSummaryFun<TData> fun)
{
    m_summary_fun = fun;
}

template<typename TData>
inline void LFMCMC<TData>::set_kernel_fun(LFMCMCKernelFun<TData> fun)
{
    m_kernel_fun = fun;
}


template<typename TData>
inline void LFMCMC<TData>::run(
    std::vector< epiworld_double > params_init_,
    size_t n_samples_,
    epiworld_double epsilon_,
    int seed
    )
{

    // Starting timing
    chrono_start();

    // Setting the baseline parameters of the model
    m_n_samples    = n_samples_;
    m_epsilon      = epsilon_;
    m_initial_params  = params_init_;
    m_n_params = params_init_.size();

    if (seed >= 0)
        this->seed(seed);

    m_current_proposed_params.resize(m_n_params);
    m_current_accepted_params.resize(m_n_params);

    if (m_simulated_data != nullptr)
        m_simulated_data->resize(m_n_samples);

    m_current_accepted_params = m_initial_params;
    m_current_proposed_params  = m_initial_params;

    // Computing the baseline sufficient statistics
    m_summary_fun(m_observed_stats, m_observed_data, this);
    m_n_stats = m_observed_stats.size();

    // Reserving size
    m_current_proposed_stats.resize(m_n_stats);
    m_current_accepted_stats.resize(m_n_stats);
    m_all_sample_drawn_prob.resize(m_n_samples);
    m_all_sample_acceptance.resize(m_n_samples, false);
    m_all_sample_params.resize(m_n_samples * m_n_params);
    m_all_sample_stats.resize(m_n_samples * m_n_stats);
    m_all_sample_kernel_scores.resize(m_n_samples);

    m_all_accepted_params.resize(m_n_samples * m_n_params);
    m_all_accepted_stats.resize(m_n_samples * m_n_stats);
    m_all_accepted_kernel_scores.resize(m_n_samples);

    TData data_i = m_simulation_fun(m_initial_params, this);

    m_summary_fun(m_current_proposed_stats, data_i, this);
    m_all_accepted_kernel_scores[0u] = m_kernel_fun(
        m_current_proposed_stats, m_observed_stats, m_epsilon, this
        );

    // Recording statistics
    for (size_t i = 0u; i < m_n_stats; ++i)
        m_all_sample_stats[i] = m_current_proposed_stats[i];
    
    m_current_accepted_stats = m_current_proposed_stats;

    for (size_t k = 0u; k < m_n_params; ++k)
        m_all_accepted_params[k] = m_initial_params[k];
    
    for (size_t k = 0u; k < m_n_params; ++k)
        m_all_sample_params[k] = m_initial_params[k];
   
    // Init progress bar
    progress_bar = Progress(m_n_samples, 80);
    if (verbose) { 
        progress_bar.next(); 
    }

    // Run LFMCMC
    for (size_t i = 1u; i < m_n_samples; ++i)
    {
        // Step 1: Generate a proposal and store it in m_current_proposed_params
        m_proposal_fun(m_current_proposed_params, m_current_accepted_params, this);

        // Step 2: Using m_current_proposed_params, simulate data
        TData data_i = m_simulation_fun(m_current_proposed_params, this);

        // Are we storing the data?
        if (m_simulated_data != nullptr)
            m_simulated_data->operator[](i) = data_i;

        // Step 3: Generate the summary statistics of the data
        m_summary_fun(m_current_proposed_stats, data_i, this);

        // Step 4: Compute the hastings ratio using the kernel function
        epiworld_double hr = m_kernel_fun(
            m_current_proposed_stats, m_observed_stats, m_epsilon, this
            );

        m_all_sample_kernel_scores[i] = hr;

        // Storing data
        for (size_t k = 0u; k < m_n_params; ++k)
            m_all_sample_params[i * m_n_params + k] = m_current_proposed_params[k];

        for (size_t k = 0u; k < m_n_stats; ++k)
            m_all_sample_stats[i * m_n_stats + k] = m_current_proposed_stats[k];
        
        // Running Hastings ratio
        epiworld_double r = runif();
        m_all_sample_drawn_prob[i] = r;

        // Step 5: Update if likely
        if (r < std::min(static_cast<epiworld_double>(1.0), hr / m_all_accepted_kernel_scores[i - 1u]))
        {
            m_all_accepted_kernel_scores[i] = hr;
            m_all_sample_acceptance[i]     = true;
            
            for (size_t k = 0u; k < m_n_stats; ++k)
                m_all_accepted_stats[i * m_n_stats + k] =
                    m_current_proposed_stats[k];

            m_current_accepted_params = m_current_proposed_params;
            m_current_accepted_stats = m_current_proposed_stats;
        } else
        {

            for (size_t k = 0u; k < m_n_stats; ++k)
                m_all_accepted_stats[i * m_n_stats + k] =
                    m_all_accepted_stats[(i - 1) * m_n_stats + k];

            m_all_accepted_kernel_scores[i] = m_all_accepted_kernel_scores[i - 1u];
        }
            

        for (size_t k = 0u; k < m_n_params; ++k)
            m_all_accepted_params[i * m_n_params + k] = m_current_accepted_params[k];

        if (verbose) { 
            progress_bar.next(); 
        }
    }

    // End timing
    chrono_end();

}


template<typename TData>
inline epiworld_double LFMCMC<TData>::runif()
{
    return runifd->operator()(*m_engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::runif(
    epiworld_double lb,
    epiworld_double ub
)
{
    return runifd->operator()(*m_engine) * (ub - lb) + lb;
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rnorm()
{
    return rnormd->operator()(*m_engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rnorm(
    epiworld_double mean,
    epiworld_double sd
    )
{
    return (rnormd->operator()(*m_engine)) * sd + mean;
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rgamma()
{
    return rgammad->operator()(*m_engine);
}

template<typename TData>
inline epiworld_double LFMCMC<TData>::rgamma(
    epiworld_double alpha,
    epiworld_double beta
    )
{

    auto old_param = rgammad->param();

    rgammad->param(std::gamma_distribution<>::param_type(alpha, beta));

    epiworld_double ans = rgammad->operator()(*m_engine);

    rgammad->param(old_param);

    return ans;

}

template<typename TData>
inline void LFMCMC<TData>::seed(epiworld_fast_uint s) {

    this->m_engine->seed(s);

}

template<typename TData>
inline void LFMCMC<TData>::set_rand_engine(std::shared_ptr< std::mt19937 > & eng)
{
    m_engine = eng;
}

template<typename TData>
inline void LFMCMC<TData>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::make_shared<std::gamma_distribution<>>(alpha,beta);
}

template<typename TData>
inline std::shared_ptr< std::mt19937 > & LFMCMC<TData>::get_rand_endgine()
{
    return m_engine;
}

// Step 1: Simulate data

// Step 2: Compute the sufficient statistics

// Step 3: Compute the hastings-ratio

// Step 4: Accept/reject, and go back to step 1

#define DURCAST(tunit,txtunit) {\
        elapsed       = std::chrono::duration_cast<std::chrono:: tunit>(\
            m_end_time - m_start_time).count(); \
        abbr_unit     = txtunit;}

template<typename TData>
inline void LFMCMC<TData>::get_elapsed_time(
    std::string unit,
    epiworld_double * last_elapsed,
    std::string * unit_abbr,
    bool print
) const {

    // Preparing the result
    epiworld_double elapsed;
    std::string abbr_unit;

    // Figuring out the length
    if (unit == "auto")
    {

        size_t tlength = std::to_string(
            static_cast<int>(floor(m_elapsed_time.count()))
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
    m_start_time = std::chrono::steady_clock::now();
}

template<typename TData>
inline void LFMCMC<TData>::chrono_end() {
    m_end_time = std::chrono::steady_clock::now();
    m_elapsed_time += (m_end_time - m_start_time);
}

template<typename TData>
inline void LFMCMC<TData>::set_params_names(std::vector< std::string > names)
{

    if (names.size() != m_n_params)
        throw std::length_error("The number of names to add differs from the number of parameters in the model.");

    m_param_names = names;

}
template<typename TData>
inline void LFMCMC<TData>::set_stats_names(std::vector< std::string > names)
{

    if (names.size() != m_n_stats)
        throw std::length_error("The number of names to add differs from the number of statistics in the model.");

    m_stat_names = names;

}

template<typename TData>
inline std::vector< epiworld_double > LFMCMC<TData>::get_mean_params()
{
    std::vector< epiworld_double > res(this->m_n_params, 0.0);
    
    for (size_t k = 0u; k < m_n_params; ++k)
    {
        for (size_t i = 0u; i < m_n_samples; ++i)
            res[k] += (this->m_all_accepted_params[k + m_n_params * i])/
                static_cast< epiworld_double >(m_n_samples);
    }

    return res;

}

template<typename TData>
inline std::vector< epiworld_double > LFMCMC<TData>::get_mean_stats()
{
    std::vector< epiworld_double > res(this->m_n_stats, 0.0);
    
    for (size_t k = 0u; k < m_n_stats; ++k)
    {
        for (size_t i = 0u; i < m_n_samples; ++i)
            res[k] += (this->m_all_accepted_stats[k + m_n_stats * i])/
                static_cast< epiworld_double >(m_n_samples);
    }

    return res;

}

template<typename TData>
inline LFMCMC<TData> & LFMCMC<TData>::verbose_off()
{
    verbose = false;
    return *this;
}

template<typename TData>
inline LFMCMC<TData> & LFMCMC<TData>::verbose_on()
{
    verbose = true;
    return *this;
}

#endif