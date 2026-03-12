#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <string_view>
#include <map>
#include <unordered_map>
#include <chrono>
#include <climits>
#include <cstdint>
#include <algorithm>
#include <regex>
#include <sstream>
#include <iomanip>
#include <set>
#include <type_traits>
#include <cassert>
#ifdef EPI_DEBUG_VIRUS
#include <atomic>
#endif

#ifndef EPIWORLD_HPP
#define EPIWORLD_HPP

/* Versioning */
#define EPIWORLD_VERSION_MAJOR 0
#define EPIWORLD_VERSION_MINOR 13
#define EPIWORLD_VERSION_PATCH 0

static const int epiworld_version_major = EPIWORLD_VERSION_MAJOR;
static const int epiworld_version_minor = EPIWORLD_VERSION_MINOR;
static const int epiworld_version_patch = EPIWORLD_VERSION_PATCH;

namespace epiworld {

    #include "config.hpp"
    #include "epiworld-macros.hpp"

    #include "misc.hpp"
    #include "progress.hpp"

    #include "rng-utils.hpp"
    #include "modeldiagram-bones.hpp"
    #include "modeldiagram-meat.hpp"

    #include "math/distributions.hpp"

    #include "math/lfmcmc.hpp"

    #include "userdata-bones.hpp"
    #include "userdata-meat.hpp"

    #include "seq_processing.hpp"

    #include "hospitalizationstracker-bones.hpp"
    #include "hospitalizationstracker-meat.hpp"

    #include "database-bones.hpp"
    #include "database-meat.hpp"
    #include "adjlist-bones.hpp"
    #include "adjlist-meat.hpp"

    #include "randgraph.hpp"

    #include "queue-bones.hpp"

    #include "globalevent-bones.hpp"
    #include "globalevent-meat.hpp"

    #include "model-bones.hpp"
    #include "model-meat.hpp"

    #include "viruses-bones.hpp"

    #include "virus-bones.hpp"
    #include "virus-distribute-meat.hpp"
    #include "virus-meat.hpp"
    
    #include "tools-bones.hpp"

    #include "tool-bones.hpp"
    #include "tool-distribute-meat.hpp"
    #include "tool-meat.hpp"

    #include "entity-bones.hpp"
    #include "entity-distribute-meat.hpp"
    #include "entity-meat.hpp"
    
    #include "agent-meat-virus-sampling.hpp"
    #include "agent-meat-state.hpp"
    #include "agent-bones.hpp"
    #include "agent-meat.hpp"

    #include "agentssample-bones.hpp"

    #include "contacttracing-bones.hpp"
    
    #include "tools/vaccine.hpp"
    #include "models/models.hpp"

}

// ---------------------------------------------------------------------------
// Specializations of std::generate_canonical for epiworld::epi_xoshiro256ss.
//
// This is the agnostic performance fix: by specializing generate_canonical for
// our custom engine, ALL std distributions (normal, gamma, binomial, etc.)
// that call generate_canonical internally will automatically bypass the
// long-double arithmetic that causes severe slowdowns on platforms where
// long double is software-emulated 128-bit (e.g., Clang+libstdc++ on ARM).
//
// The specialisation converts the 64-bit engine output to a floating-point
// value in [0, 1) using only the precision of the target type:
//   - double (53 mantissa bits): top 53 bits of the 64-bit word, × 2^{-53}
//   - float  (24 mantissa bits): top 24 bits,                      × 2^{-24}
//
// The result is strictly less than 1.0 by construction (max value is
// (2^N - 1) / 2^N), so no nextafter clamp is required.
// ---------------------------------------------------------------------------
namespace std {
    template<>
    inline double generate_canonical<
        double,
        numeric_limits<double>::digits,
        epiworld::epi_xoshiro256ss
    >(epiworld::epi_xoshiro256ss& eng) {
        return static_cast<double>(eng() >> 11) * 0x1.0p-53;
    }

    template<>
    inline float generate_canonical<
        float,
        numeric_limits<float>::digits,
        epiworld::epi_xoshiro256ss
    >(epiworld::epi_xoshiro256ss& eng) {
        return static_cast<float>(eng() >> 40) * 0x1.0p-24f;
    }
}

#endif
