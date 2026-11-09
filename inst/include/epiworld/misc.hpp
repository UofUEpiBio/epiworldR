#ifndef EPIWORLD_MISC_HPP
#define EPIWORLD_MISC_HPP

class Model;

class Agent;

// Relevant for anything using vecHasher function ------------------------------
/**
 * @brief Vector hasher
 * @tparam T
 */
template <typename T>
struct vecHasher {
    std::size_t operator()(std::vector< T > const&  dat) const noexcept {

        std::hash< T > hasher;
        std::size_t hash = hasher(dat[0u]);

        // ^ makes bitwise XOR
        // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
        // << is a shift operator, something like lhs * 2^(rhs)
        if (dat.size() > 1u)
            for (epiworld_fast_uint i = 1u; i < dat.size(); ++i)
                hash ^= hasher(dat[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);

        return hash;

    }
};

template<typename Ta = epiworld_double, typename Tb = epiworld_fast_uint>
using MapVec_type = std::unordered_map< std::vector< Ta >, Tb, vecHasher<Ta>>;

/**
 * @brief Check whether `a` is included in `b`
 *
 * @tparam Ta Type of `a`. Could be int, epiworld_double, etc.
 * @param a Scalar of class `Ta`.
 * @param b Vector `std::vector` of class `Ta`.
 * @return `true` if `a in b`, and `false` otherwise.
 */
template<typename Ta>
inline bool IN(const Ta & a, const std::vector< Ta > & b) noexcept
{
    for (const auto & i : b)
        if (a == i)
            return true;

    return false;
}

/**
 * @brief Read parameters from a yaml file
 *
 * @details
 * The file should have the following structure:
 * ```yaml
 * # Comment
 * [name of parameter 1]: [value in T]
 * [name of parameter 2]: [value in T]
 * ...
 * ```
 *
 * @tparam T Type of the parameter
 * @param fn Path to the file containing the parameters
 * @return std::map<std::string, T>
 */
template <typename T>
inline std::map< std::string, T > read_yaml(std::string fn)
{

    std::ifstream paramsfile(fn);

    if (!paramsfile)
        throw std::logic_error("The file " + fn + " was not found.");

    std::regex pattern("^([^:]+)\\s*[:]\\s*([-]?[0-9]+|[-]?[0-9]*\\.[0-9]+)?\\s*$");

    std::string line;
    std::smatch match;
    auto empty = std::sregex_iterator();

    // Making room
    std::map<std::string, T> parameters;

    while (std::getline(paramsfile, line))
    {

        // Is it a comment or an empty line?
        if (std::regex_match(line, std::regex("^([*].+|//.+|#.+|\\s*)$")))
            continue;

        // Finding the pattern, if it doesn't match, then error
        std::regex_match(line, match, pattern);

        if (match.empty())
            throw std::logic_error("Line has invalid format:\n" + line);

        // Capturing the number
        std::string anumber = match[2u].str() + match[3u].str();
        T tmp_num = static_cast<T>(
            std::strtod(anumber.c_str(), nullptr)
            );

        std::string pname = std::regex_replace(
            match[1u].str(),
            std::regex("^\\s+|\\s+$"),
            "");

        // Adding the parameter to the map
        parameters[pname] = tmp_num;

    }

    return parameters;

}



/**
 * @brief Cast a `Model` pointer to a specific model type
 */
#ifndef EPI_ASSUME
    #if defined(EPI_DEBUG)
        #include <cassert>
        #define EPI_ASSUME(cond) assert(cond)
    #elif defined(__clang__)
        #define EPI_ASSUME(cond) __builtin_assume(cond)
    #elif defined(_MSC_VER)
        #define EPI_ASSUME(cond) __assume(cond)
    #elif defined(__GNUC__)
        #define EPI_ASSUME(cond) do { if (!(cond)) __builtin_unreachable(); } while (0)
    #else
        #define EPI_ASSUME(cond) ((void)0)
    #endif
#endif

/**
 * @brief Index of the lowest set bit of a non-zero 64-bit word.
 *
 * @details Used to walk bitsets of agents (the queue's queued agents, the
 * carriers that push, the agents to update after a push) in ascending id
 * order: take the lowest set bit, visit that agent, clear the bit, repeat. It
 * costs one instruction per agent visited, and a whole word of 64 absent
 * agents is skipped at once.
 *
 * C++17 has no `std::countr_zero` (that is C++20). GCC and Clang get the
 * builtin; other compilers use a de Bruijn multiply, which needs no platform
 * header. `x` must not be zero.
 */
inline unsigned int epi_ctz64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return static_cast< unsigned int >(__builtin_ctzll(x));
#else
    static const unsigned int table[64] = {
         0,  1, 48,  2, 57, 49, 28,  3, 61, 58, 50, 42, 38, 29, 17,  4,
        62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12,  5,
        63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
        46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19,  9, 13,  8,  7,  6
    };
    return table[((x & (~x + 1u)) * 0x03f79d71b4cb0a89ull) >> 58];
#endif
}

template<class To>
inline To* model_cast(Model* m) {
#ifdef EPI_DEBUG
    auto* out = dynamic_cast<To*>(m);
    assert(out != nullptr);
    return out;
#else
    auto* out = static_cast<To*>(m);
    EPI_ASSUME(out != nullptr);
    return out;
#endif
}

#endif
