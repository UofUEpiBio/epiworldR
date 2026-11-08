#ifndef EPIWORLD_MODEL_RAND_MEAT_HPP
#define EPIWORLD_MODEL_RAND_MEAT_HPP

#include "rng-utils.hpp"
#include "model-bones.hpp"

inline void Model::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::gamma_distribution<>(alpha,beta);
}

inline void Model::set_rand_norm(epiworld_double mean, epiworld_double sd)
{
    rnormd  = std::normal_distribution<>(mean, sd);
}

inline void Model::set_rand_unif(epiworld_double a, epiworld_double b)
{
    runifd_a = a;
    runifd_b = b;
}

inline void Model::set_rand_lognormal(epiworld_double mean, epiworld_double shape)
{
    rlognormald  = std::lognormal_distribution<>(mean, shape);
}

inline void Model::set_rand_exp(epiworld_double lambda)
{
    rexpd  = std::exponential_distribution<>(lambda);
}

inline void Model::set_rand_binom(int n, epiworld_double p)
{
    rbinomd  = std::binomial_distribution<>(n, p);
    rbinomd_n = n;
#ifdef EPI_FAST_BINOM
    rbinomd_fast_lambda = static_cast<epiworld_double>(n) * p;
    // A Poisson needs a positive mean; with n = 0 or p = 0 the binomial
    // gives 0 anyway
    rbinomd_use_poisson =
        (rbinomd_fast_lambda > 0.0) &&
        (p <= static_cast<epiworld_double>(0.01)) &&
        (static_cast<epiworld_double>(n) * p * p <= static_cast<epiworld_double>(0.1));
#else
    rbinomd_fast_lambda = 0.0;
    rbinomd_use_poisson = false;
#endif
}

inline void Model::set_rand_nbinom(int n, epiworld_double p)
{
    rnbinomd  = std::negative_binomial_distribution<>(n, p);
}

inline void Model::set_rand_geom(epiworld_double p)
{
    rgeomd  = std::geometric_distribution<>(p);
}

inline void Model::set_rand_poiss(epiworld_double lambda)
{
    rpoissd  = std::poisson_distribution<>(lambda);
}

inline std::shared_ptr< epi_xoshiro256ss > & Model::get_rand_endgine()
{
    return engine;
}

inline void Model::set_rand_engine(std::shared_ptr< epi_xoshiro256ss > & eng)
{
    engine = eng;
}

inline epiworld_double Model::runif() {
    // CHECK_INIT()
    epiworld_double res = runif_epi(*engine);
    return res * (runifd_b - runifd_a) + runifd_a;
}

inline int Model::runif_int(int a, int b) {
    // CHECK_INIT()
    return a + runif_index(static_cast<uint32_t>(b - a + 1));
}

inline uint32_t Model::runif_index(uint32_t n) {

    // If asking for a 0-range, return 0 to prevent division by zero in the bias check
    if (n == 0) return 0;

    // Grab 32 perfectly uniform random bits directly from xoshiro256ss
    uint32_t x = static_cast<uint32_t>((*engine)());
    
    // Multiply by the bound N to get a 64-bit result
    uint64_t m = static_cast<uint64_t>(x) * static_cast<uint64_t>(n);
    
    // The fractional part (lower 32 bits)
    uint32_t l = static_cast<uint32_t>(m);
    
    // Unbiased rejection sampling for the rare edge case
    if (l < n) {
        uint32_t t = -n % n; // Two's complement trick to get (2^32 - n) % n
        while (l < t) {
            x = static_cast<uint32_t>((*engine)());
            m = static_cast<uint64_t>(x) * static_cast<uint64_t>(n);
            l = static_cast<uint32_t>(m);
        }
    }
    
    // The integer part (upper 32 bits) is our unbiased scaled random index
    return static_cast<uint32_t>(m >> 32);
}

inline epiworld_double Model::runif(epiworld_double a, epiworld_double b) {
    // CHECK_INIT()
    return runif_epi(*engine) * (b - a) + a;
}

inline epiworld_double Model::rnorm() {
    // CHECK_INIT()
    return rnormd(*engine);
}

inline epiworld_double Model::rnorm(epiworld_double mean, epiworld_double sd) {
    // CHECK_INIT()
    return rnormd(*engine) * sd + mean;
}

inline epiworld_double Model::rgamma() {
    return rgammad(*engine);
}

inline epiworld_double Model::rgamma(epiworld_double alpha, epiworld_double beta) {

    return rgammad(
        *engine,
        std::gamma_distribution<>::param_type(alpha, beta)
    );

}

inline epiworld_double Model::rexp() {
    return rexpd(*engine);
}

inline epiworld_double Model::rexp(epiworld_double lambda) {

    return rexpd(
        *engine,
        std::exponential_distribution<>::param_type(lambda)
    );

}

inline epiworld_double Model::rlognormal() {
    return rlognormald(*engine);
}

inline epiworld_double Model::rlognormal(epiworld_double mean, epiworld_double shape) {

    return rlognormald(
        *engine,
        std::lognormal_distribution<>::param_type(mean, shape)
    );
}

inline int Model::rbinom() {
#ifdef EPI_FAST_BINOM
    if (rbinomd_use_poisson)
    {
        int res = rpoiss(rbinomd_fast_lambda);
        return std::min(res, rbinomd_n);
    }
#endif
    return rbinomd(*engine);
}

inline int Model::rbinom(int n, epiworld_double p) {

    if (n == 0 || p == 0.0)
        return 0;

#ifdef EPI_FAST_BINOM
    if (
        (p <= static_cast<epiworld_double>(0.01)) &&
        (static_cast<epiworld_double>(n) * p * p <= static_cast<epiworld_double>(0.1))
    )
    {
        int res = rpoiss(static_cast<epiworld_double>(n) * p);
        return std::min(res, n);
    }
#endif

    return rbinomd(
        *engine,
        std::binomial_distribution<>::param_type(n, p)
    );

}

inline int Model::rnbinom() {
    return rnbinomd(*engine);
}

inline int Model::rnbinom(int n, epiworld_double p) {

    return rnbinomd(
        *engine,
        std::negative_binomial_distribution<>::param_type(n, p)
    );
}

inline int Model::rgeom() {
    return rgeomd(*engine);
}

inline int Model::rgeom(epiworld_double p) {

    return rgeomd(
        *engine,
        std::geometric_distribution<>::param_type(p)
    );

}

inline int Model::rpoiss() {
    return rpoissd(*engine);
}

inline int Model::rpoiss(epiworld_double lambda) {

    return rpoissd(
        *engine,
        std::poisson_distribution<>::param_type(lambda)
    );

}

inline size_t Model::sample_from_probs(size_t n) {

    epiworld_double p_total = runif();
    size_t ans;
    for (ans = 0u; ans < n; ++ans)
    {
        if (p_total < array_double_tmp[ans])
            break;
        array_double_tmp[ans + 1] += array_double_tmp[ans];
    }
    return ans;

}

inline void Model::reserve_tmp(size_t n) {

    if (n <= array_virus_tmp.size())
        return;

    // Doubling keeps the growth amortized
    n = std::max(n, 2u * array_virus_tmp.size());
    array_virus_tmp.resize(n);
    array_layer_tmp.resize(n);
    array_double_tmp.resize(2u * n);

}

inline void Model::seed(size_t s) {

    this->engine->seed(s);

    // Some distributions keep state between calls (e.g., the normal draws
    // made in pairs, also inside the binomial); without a reset, a run would
    // depend on the runs before it, not only on its seed
    rnormd.reset();
    rgammad.reset();
    rlognormald.reset();
    rexpd.reset();
    rbinomd.reset();
    rnbinomd.reset();
    rgeomd.reset();
    rpoissd.reset();

}

/**
 * @brief Conditional Weighted Sampling
 *
 * @details
 * The sampling function will draw one of `{-1, 0,...,probs.size() - 1}` in a
 * weighted fashion. The probabilities are drawn given that either one or none
 * of the cases is drawn; in the latter returns -1.
 *
 * @param probs Vector of probabilities.
 * @param m A `Model`. This is used to draw random uniform numbers.
 * @return int If -1 then it means that none got sampled, otherwise the index
 * of the entry that got drawn.
 */
template<typename TDbl = epiworld_double>
inline int roulette(
    const std::vector< TDbl > & probs,
    Model * m
    )
{

    // Step 1: a probability of 1 is a certain event; if there are any, one of
    // them is drawn uniformly. (This used to test `> 1 - 1e-100`, which is
    // `> 1` in floating point, so p == 1 fell through to 0/0 below and the
    // last entry won regardless.)
    std::vector< int > certain_infection;
    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
        if (probs[p] >= 1.0)
            certain_infection.push_back(static_cast<int>(p));

    TDbl r = static_cast<TDbl>(m->runif());
    if (certain_infection.size() > 0)
        return certain_infection[std::floor(r * certain_infection.size())];

    // Step 2: the odds r_p = p / (1 - p) and their sum R. Then
    // P(none) = 1 / (1 + R) and P(p) = r_p / (1 + R), the same as normalizing
    // by prod(1 - p), which underflows to 0 with many entries (0/0 weights).
    std::vector< TDbl > odds(probs.size());
    TDbl odds_total = 0.0;
    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        odds[p] = probs[p] / (1.0 - probs[p]);
        odds_total += odds[p];
    }

    // Step 3: Roulette
    TDbl cumsum = 1.0 / (1.0 + odds_total);
    if (r < cumsum)
        return -1;

    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        // If it yield here, then bingo, the individual will acquire the disease
        cumsum += odds[p] / (1.0 + odds_total);
        if (r < cumsum)
            return static_cast<int>(p);

    }

    #ifdef EPI_DEBUG
    printf_epiworld("[epi-debug] roulette::cumsum = %.4f\n", cumsum);
    #endif

    return static_cast<int>(probs.size() - 1u);

}

inline int roulette(std::vector< double > & probs, Model * m)
{
    return roulette<double>(probs, m);
}

inline int roulette(std::vector< float > & probs, Model * m)
{
    return roulette<float>(probs, m);
}


inline int roulette(
    epiworld_fast_uint nelements,
    Model * m
    )
{

    if ((nelements * 2) > m->array_double_tmp.size())
    {
        throw std::logic_error(
            "Trying to sample from more data than there is in roulette!" +
            std::to_string(nelements) + " vs " +
            std::to_string(m->array_double_tmp.size())
            );
    }

    auto & tmp = m->array_double_tmp;

    // Step 1: certain events (p >= 1; see the vector version), drawn uniformly
    epiworld_fast_uint ncertain = 0u;
    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
        if (tmp[p] >= 1.0)
            tmp[nelements + ncertain++] = p;

    epiworld_double r = m->runif();
    if (ncertain > 0u)
        return tmp[nelements + std::floor(ncertain * r)];

    // Step 2: the odds and their sum, normalized by 1 + sum (see the vector
    // version; prod(1 - p) underflows with many exposures)
    epiworld_double odds_total = 0.0;
    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
    {
        tmp[nelements + p] = tmp[p] / (1.0 - tmp[p]);
        odds_total += tmp[nelements + p];
    }

    // Step 3: Roulette
    epiworld_double cumsum = 1.0 / (1.0 + odds_total);
    if (r < cumsum)
        return -1;

    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
    {
        // If it yield here, then bingo, the individual will acquire the disease
        cumsum += tmp[nelements + p] / (1.0 + odds_total);
        if (r < cumsum)
            return static_cast<int>(p);

    }

    return static_cast<int>(nelements - 1u);

}

#endif
