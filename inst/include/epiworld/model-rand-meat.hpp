#ifndef EPIWORLD_MODEL_RAND_MEAT_HPP
#define EPIWORLD_MODEL_RAND_MEAT_HPP

#include "rng-utils.hpp"
#include "model-bones.hpp"

template<typename TSeq>
inline void Model<TSeq>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::gamma_distribution<>(alpha,beta);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_norm(epiworld_double mean, epiworld_double sd)
{
    rnormd  = std::normal_distribution<>(mean, sd);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_unif(epiworld_double a, epiworld_double b)
{
    runifd_a = a;
    runifd_b = b;
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_lognormal(epiworld_double mean, epiworld_double shape)
{
    rlognormald  = std::lognormal_distribution<>(mean, shape);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_exp(epiworld_double lambda)
{
    rexpd  = std::exponential_distribution<>(lambda);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_binom(int n, epiworld_double p)
{
    rbinomd  = std::binomial_distribution<>(n, p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_nbinom(int n, epiworld_double p)
{
    rnbinomd  = std::negative_binomial_distribution<>(n, p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_geom(epiworld_double p)
{
    rgeomd  = std::geometric_distribution<>(p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_poiss(epiworld_double lambda)
{
    rpoissd  = std::poisson_distribution<>(lambda);
}

template<typename TSeq>
inline std::shared_ptr< epi_xoshiro256ss > & Model<TSeq>::get_rand_endgine()
{
    return engine;
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_engine(std::shared_ptr< epi_xoshiro256ss > & eng)
{
    engine = eng;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif() {
    // CHECK_INIT()
    epiworld_double res = runif_epi(*engine);
    return res * (runifd_b - runifd_a) + runifd_a;
}

template<typename TSeq>
inline int Model<TSeq>::runif_int(int a, int b) {
    // CHECK_INIT()
    auto res =
        static_cast<int>(std::floor(runif_epi(*engine) * (b - a + 1))) + 
        a;

    // Checking it is within the bounds
    if (res < a)
        res = a;
    else if (res > b)
        res = b;

    return res;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif(epiworld_double a, epiworld_double b) {
    // CHECK_INIT()
    return runif_epi(*engine) * (b - a) + a;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm() {
    // CHECK_INIT()
    return rnormd(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm(epiworld_double mean, epiworld_double sd) {
    // CHECK_INIT()
    return rnormd(*engine) * sd + mean;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma() {
    return rgammad(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma(epiworld_double alpha, epiworld_double beta) {

    return rgammad(
        *engine,
        std::gamma_distribution<>::param_type(alpha, beta)
    );

}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp() {
    return rexpd(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp(epiworld_double lambda) {

    return rexpd(
        *engine,
        std::exponential_distribution<>::param_type(lambda)
    );

}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal() {
    return rlognormald(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal(epiworld_double mean, epiworld_double shape) {

    return rlognormald(
        *engine,
        std::lognormal_distribution<>::param_type(mean, shape)
    );
}

template<typename TSeq>
inline int Model<TSeq>::rbinom() {
    return rbinomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rbinom(int n, epiworld_double p) {

    if (n == 0 || p == 0.0)
        return 0;

    return rbinomd(
        *engine,
        std::binomial_distribution<>::param_type(n, p)
    );

}

template<typename TSeq>
inline int Model<TSeq>::rnbinom() {
    return rnbinomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rnbinom(int n, epiworld_double p) {

    return rnbinomd(
        *engine,
        std::negative_binomial_distribution<>::param_type(n, p)
    );
}

template<typename TSeq>
inline int Model<TSeq>::rgeom() {
    return rgeomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rgeom(epiworld_double p) {

    return rgeomd(
        *engine,
        std::geometric_distribution<>::param_type(p)
    );

}

template<typename TSeq>
inline int Model<TSeq>::rpoiss() {
    return rpoissd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rpoiss(epiworld_double lambda) {

    return rpoissd(
        *engine,
        std::poisson_distribution<>::param_type(lambda)
    );

}

template<typename TSeq>
inline size_t Model<TSeq>::sample_from_probs(size_t n) {

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

template<typename TSeq>
inline void Model<TSeq>::seed(size_t s) {
    this->engine->seed(s);
}

#endif