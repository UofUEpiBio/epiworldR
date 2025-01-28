#ifndef LFMCMC_MEAT_PRINT_HPP
#define LFMCMC_MEAT_PRINT_HPP

template<typename TData>
inline void LFMCMC<TData>::print(size_t burnin) const
{

    // For each statistic or parameter in the model, we print three values: 
    // - mean, the 2.5% quantile, and the 97.5% quantile
    std::vector< epiworld_double > summ_params(m_n_params * 3, 0.0);
    std::vector< epiworld_double > summ_stats(m_n_stats * 3, 0.0);

    // Compute the number of samples to use based on burnin rate
    size_t n_samples_print = m_n_samples;
    if (burnin > 0)
    {
        if (burnin >= m_n_samples)
            throw std::length_error(
                "The burnin is greater than or equal to the number of samples."
                );

        n_samples_print = m_n_samples - burnin;

    }

    epiworld_double n_samples_dbl = static_cast< epiworld_double >(
        n_samples_print
        );

    // Compute parameter summary values
    for (size_t k = 0u; k < m_n_params; ++k)
    {

        // Retrieving the relevant parameter
        std::vector< epiworld_double > par_i(n_samples_print);
        for (size_t i = burnin; i < m_n_samples; ++i)
        {
            par_i[i-burnin] = m_all_accepted_params[i * m_n_params + k];
            summ_params[k * 3] += par_i[i-burnin]/n_samples_dbl;
        }

        // Computing the 95% Credible interval
        std::sort(par_i.begin(), par_i.end());

        summ_params[k * 3 + 1u] = 
            par_i[std::floor(.025 * n_samples_dbl)];
        summ_params[k * 3 + 2u] = 
            par_i[std::floor(.975 * n_samples_dbl)];

    }

    // Compute statistics summary values
    for (size_t k = 0u; k < m_n_stats; ++k)
    {

        // Retrieving the relevant parameter
        std::vector< epiworld_double > stat_k(n_samples_print);
        for (size_t i = burnin; i < m_n_samples; ++i)
        {
            stat_k[i-burnin] = m_all_accepted_stats[i * m_n_stats + k];
            summ_stats[k * 3] += stat_k[i-burnin]/n_samples_dbl;
        }

        // Computing the 95% Credible interval
        std::sort(stat_k.begin(), stat_k.end());
        summ_stats[k * 3 + 1u] = 
            stat_k[std::floor(.025 * n_samples_dbl)];
        summ_stats[k * 3 + 2u] = 
            stat_k[std::floor(.975 * n_samples_dbl)];

    }

    printf_epiworld("___________________________________________\n\n");
    printf_epiworld("LIKELIHOOD-FREE MARKOV CHAIN MONTE CARLO\n\n");

    printf_epiworld("N Samples (total) : %zu\n", m_n_samples);
    printf_epiworld("N Samples (after burn-in period) : %zu\n", m_n_samples - burnin);

    std::string abbr;
    epiworld_double elapsed;
    get_elapsed_time("auto", &elapsed, &abbr, false);
    printf_epiworld("Elapsed t : %.2f%s\n\n", elapsed, abbr.c_str());
    
    ////////////////////////////////////////////////////////////////////////////
    // PARAMETERS
    ////////////////////////////////////////////////////////////////////////////
    printf_epiworld("Parameters:\n");

    // Figuring out format
    std::string fmt_params;
    
    int nchar_par_num = 0;
    for (auto & n : summ_params)
    {
        
        int tmp_nchar;
        
        if (std::abs(n) < 1) {
            // std::log10(<1) will return negative number
            // std::log10(0) will return -inf and throw a runtime error
            tmp_nchar = 0;
        } else {
            tmp_nchar = std::floor(std::log10(std::abs(n)));
        }

        if (nchar_par_num < tmp_nchar)
            nchar_par_num = tmp_nchar;
    }
    nchar_par_num += 5; // 1 for neg padd, 2 for decimals, 1 the decimal point, and one b/c log(<10) < 1.
    std::string charlen = std::to_string(nchar_par_num);

    if (m_param_names.size() != 0u)
    {
        int nchar_par = 0;
        for (auto & n : m_param_names)
        {
            int tmp_nchar = n.length();
            if (nchar_par < tmp_nchar)
                nchar_par = tmp_nchar;
        }

        fmt_params = std::string("  -%-") +
            std::to_string(nchar_par) +
            std::string("s : % ") + charlen  + 
            std::string(".2f [% ") + charlen + 
            std::string(".2f, % ") + charlen +
            std::string(".2f] (initial : % ") +
            charlen + std::string(".2f)\n");

        for (size_t k = 0u; k < m_n_params; ++k)
        {
            printf_epiworld(
                fmt_params.c_str(),
                m_param_names[k].c_str(),
                summ_params[k * 3],
                summ_params[k * 3 + 1u],
                summ_params[k * 3 + 2u],
                m_initial_params[k]
                );
        }

        
    } else {

        fmt_params = std::string("  [%-2ld]: % ") + charlen + 
            std::string(".2f [% ") + charlen +
            std::string(".2f, % ") + charlen + 
            std::string(".2f] (initial : % ") + charlen +
            std::string(".2f)\n");

        for (size_t k = 0u; k < m_n_params; ++k)
        {
            
            printf_epiworld(
                fmt_params.c_str(),
                k,
                summ_params[k * 3],
                summ_params[k * 3 + 1u],
                summ_params[k * 3 + 2u],
                m_initial_params[k]
                );
        }

    }    

    ////////////////////////////////////////////////////////////////////////////
    // Statistics
    ////////////////////////////////////////////////////////////////////////////
    printf_epiworld("\nStatistics:\n");
    int nchar = 0;
    for (auto & s : summ_stats)
    {
        int tmp_nchar;
        if (std::abs(s) < 1) {
            // std::log10(<1) will return negative number
            // std::log10(0) will return -inf and throw a runtime error
            tmp_nchar = 0;
        } else {
            tmp_nchar = std::floor(std::log10(std::abs(s)));
        }
    
        if (nchar < tmp_nchar)
            nchar = tmp_nchar;
    }

    nchar += 5; // See above

    std::string nchar_char = std::to_string(nchar);

    // Figuring out format
    std::string fmt_stats;
    if (m_stat_names.size() != 0u)
    {
        int nchar_stats = 0;
        for (auto & n : m_stat_names)
        {
            int tmp_nchar = n.length();
            if (nchar_stats < tmp_nchar)
                nchar_stats = tmp_nchar;
        }

        fmt_stats = std::string("  -%-") +
            std::to_string(nchar_stats) +
            std::string("s : % ") + nchar_char +
            std::string(".2f [% ") + nchar_char +
            std::string(".2f, % ") + nchar_char +
            std::string(".2f] (Observed: % ") + nchar_char +
            std::string(".2f)\n");

        for (size_t k = 0u; k < m_n_stats; ++k)
        {
            printf_epiworld(
                fmt_stats.c_str(),
                m_stat_names[k].c_str(),
                summ_stats[k * 3],
                summ_stats[k * 3 + 1u],
                summ_stats[k * 3 + 2u],
                m_observed_stats[k]
                );
        }

        
    } else {

        fmt_stats = std::string("  [%-2ld] : % ") + 
            nchar_char +
            std::string(".2f [% ") + nchar_char +
            std::string(".2f, % ") + nchar_char +
            std::string(".2f] (Observed: % ") + nchar_char +
            std::string(".2f)\n");

        for (size_t k = 0u; k < m_n_stats; ++k)
        {
            printf_epiworld(
                fmt_stats.c_str(),
                k,
                summ_stats[k * 3],
                summ_stats[k * 3 + 1u],
                summ_stats[k * 3 + 2u],
                m_observed_stats[k]
                );
        }

    }

    printf_epiworld("___________________________________________\n\n");
}

#endif