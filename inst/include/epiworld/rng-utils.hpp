#ifndef EPIWORLD_RNG_UTILS_HPP
#define EPIWORLD_RNG_UTILS_HPP

#include <cstdint>
#include <limits>
#include <random>
#include "config.hpp"

/**
 * @brief Fast 64-bit PRNG based on the xoshiro256** algorithm.
 *
 * xoshiro256** has 256-bit state, 64-bit output, and passes all known
 * statistical tests. It is significantly faster than std::mt19937.
 * State is seeded via SplitMix64 from a single 64-bit value.
 *
 * Reference: Blackman & Vigna, "Scrambled Linear Pseudorandom Number
 * Generators", ACM TOMS, 2021. https://prng.di.unimi.it/
 */
class epi_xoshiro256ss {
    uint64_t s[4];

    static uint64_t rotl(const uint64_t x, int k) noexcept {
        return (x << k) | (x >> (64 - k));
    }

    static uint64_t splitmix64(uint64_t& x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

public:
    using result_type = uint64_t;

    explicit epi_xoshiro256ss(uint64_t seed_val = 0) noexcept {
        seed(seed_val);
    }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept {
        return std::numeric_limits<result_type>::max();
    }

    void seed(uint64_t seed_val) noexcept {
        s[0] = splitmix64(seed_val);
        s[1] = splitmix64(seed_val);
        s[2] = splitmix64(seed_val);
        s[3] = splitmix64(seed_val);
    }

    result_type operator()() noexcept {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }
};

/**
 * @brief Draw a uniform [0, 1) random number from an epi_xoshiro256ss engine.
 *
 * Uses the top N bits of the 64-bit output (where N =
 * std::numeric_limits<epiworld_double>::digits) to produce a value in [0, 1)
 * without any long-double arithmetic.
 *
 * @param engine An epi_xoshiro256ss engine.
 * @return epiworld_double in [0, 1).
 */
inline epiworld_double runif_epi(epi_xoshiro256ss & engine) {
    static_assert(
        std::numeric_limits<epiworld_double>::digits < 64,
        "epiworld_double must have fewer than 64 mantissa bits; "
        "the bit-extraction in runif_epi requires digits < 64 to avoid "
        "undefined behaviour in the shift and scale computation"
    );
    
    constexpr int bits  = std::numeric_limits<epiworld_double>::digits;
    constexpr int shift = 64 - bits;
    // scale = 2^{-bits}: the result is in [0, 1) by construction (no clamp needed)
    constexpr epiworld_double scale =
        epiworld_double(1) / epiworld_double(uint64_t(1) << bits);
    return static_cast<epiworld_double>(engine() >> shift) * scale;
}

/**
 * @brief Draw a uniform [0, 1) random number from a std::mt19937 engine.
 *
 * Kept for backward compatibility with code that still holds a raw mt19937.
 */
inline epiworld_double runif_mt19937(std::mt19937 & engine) {
    static constexpr epiworld_double factor =
        epiworld_double(1) / (epiworld_double(std::mt19937::max()) + epiworld_double(1));
    epiworld_double res = epiworld_double(engine()) * factor;
    if (res >= epiworld_double(1))
        res = std::nextafter(epiworld_double(1), epiworld_double(0));
    return res;
}

#endif
