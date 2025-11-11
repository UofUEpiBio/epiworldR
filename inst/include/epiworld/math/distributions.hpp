#ifndef EPIWORLD_MATH_DISTRIBUTIONS_HPP
#define EPIWORLD_MATH_DISTRIBUTIONS_HPP

// Implementing the factorial function
/**
 * @brief Compute the log of the factorial
 * 
 * @param n Number
 * 
 * @return The log of the factorial
 */
inline double log_factorial(int n)
{
    if (n == 0)
        return 0.0;
    return std::log(static_cast<double>(n)) + log_factorial(n-1);
}

/**
 * @brief Compute the Poisson probability
 * 
 * @param k Number of events
 * @param lambda Rate
 * @param max_n Maximum number of events
 * @param as_log Return the log of the probability
 * 
 * @return The Poisson probability
 */
inline double dpois(
    int k,
    double lambda,
    int max_n = 100,
    bool as_log = false
    )
{

    if (max_n < k)
        throw std::runtime_error("max_n must be greater than k");

    double res = k * std::log(lambda) - lambda - log_factorial(
        std::min(k, max_n)
        );
    
    return as_log ? res : std::exp(res);
}

/**
 * @brief Compute the probability of the generation interval
 * 
 * @details
 * If `p_0_approx` is negative, it will be computed using the Poisson
 * distribution. If `normalizing` is negative, it will be computed on the fly
 * 
 * @param g Generation interval
 * @param S Population size
 * @param p_c Probability of contact
 * @param p_i Probability of infection
 * @param p_r Probability of recovery
 * @param p_0_approx Approximation of the probability of not being infected
 * @param normalizing Normalizing constant
 * @param max_contacts Maximum number of contacts
 * @param max_days Maximum number of days
 * 
 * @return The probability of the generation interval
 * 
 */
inline double dgenint(
    int g,
    double S,
    double p_c,
    double p_i,
    double p_r,
    double & p_0_approx,
    double & normalizing,
    int max_contacts = 200,
    int max_days = 200
    ) {

    if ((g < 1) || (g > max_days))
        return 0.0;

    if (p_0_approx < 0.0)
    {

        p_0_approx = 0.0;
        for (int i = 0; i < max_contacts; ++i)
        {

            p_0_approx += std::exp(
                dpois(i, S * p_c, max_contacts, true) +
                std::log(1.0 - p_i) * static_cast<double>(i)
                );

        }
    }

    double g_dbl = static_cast<double>(g);

    if (normalizing < 0.0)
    {

        normalizing = 1.0;
        double log1_p_r = std::log(1.0 - p_r);
        double log_p_r = std::log(p_r);
        double log_p_0_approx = std::log(p_0_approx);
        for (int i = 1; i <= max_days; ++i)
        {

            double i_dbl = static_cast<double>(i);

            normalizing -= std::exp(
                log1_p_r * (i_dbl - 1.0) +
                log_p_r +
                log_p_0_approx * (i_dbl - 1.0)
                );
        }

    }


    return std::exp(
        std::log(1 - p_r) * (g_dbl)+
        std::log(p_0_approx) * (g_dbl - 1.0) +
        std::log(1.0 - p_0_approx) -
        std::log(normalizing)
        );

}

// Mean of the generation interval
/**
 * @brief Compute the mean of the generation interval
 * @param S Population size.
 * @param p_c Probability of contact.
 * @param p_i Probability of infection.
 * @param p_r Probability of recovery.
 * @param max_days Maximum number of days.
 * @param max_n Maximum number of contacts.
 * 
 * @return The mean of the generation interval
 */
inline double gen_int_mean(
    double S,
    double p_c,
    double p_i,
    double p_r,
    int max_days = 200,
    int max_n = 200
    ) {

    double mean = 0.0;
    double p_0_approx = -1.0;
    double normalizing = -1.0;
    for (int i = 1; i < max_days; ++i)
    {
        mean += 
            static_cast<double>(i) *
            dgenint(
                i, S, p_c, p_i, p_r, p_0_approx, normalizing, max_n, max_days
                );

    }

    return mean;

}

#endif