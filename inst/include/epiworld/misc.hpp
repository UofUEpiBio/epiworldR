#ifndef EPIWORLD_MISC_HPP 
#define EPIWORLD_MISC_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
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
 * @name Default sequence initializers
 * 
 * @details 
 * If the user does not provide a default sequence, this function is used when
 * a sequence needs to be initialized. Some examples: `Agent`, `Virus`, and
 * `Tool` need a default sequence.
 * 
 * @tparam TSeq 
 * @return TSeq 
 */
///@{
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline TSeq default_sequence(int seq_count);

// Making it 'static' so that we don't have problems when including the
// header. This is important during the linkage, e.g., in R.
// See https://en.cppreference.com/w/cpp/language/storage_duration#Linkage
// static int _n_sequences_created = 0;

template<>
inline bool default_sequence(int seq_count) {

    if (seq_count == 2)
        throw std::logic_error("Maximum number of sequence created.");

    return seq_count++ ? false : true;
}

template<>
inline int default_sequence(int seq_count) {
    return seq_count++;
}

template<>
inline epiworld_double default_sequence(int seq_count) {
    return static_cast<epiworld_double>(seq_count++);
}

template<>
inline std::vector<bool> default_sequence(int seq_count) {

    if (seq_count == 2)
        throw std::logic_error("Maximum number of sequence created.");

    return {seq_count++ ? false : true};
}

template<>
inline std::vector<int> default_sequence(int seq_count) {
    return {seq_count++};
}

template<>
inline std::vector<epiworld_double> default_sequence(int seq_count) {
    return {static_cast<epiworld_double>(seq_count++)};
}
///@}

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
template<typename TSeq = EPI_DEFAULT_TSEQ, typename TDbl = epiworld_double >
inline int roulette(
    const std::vector< TDbl > & probs,
    Model<TSeq> * m
    )
{

    // Step 1: Computing the prob on none 
    TDbl p_none = 1.0;
    std::vector< int > certain_infection;
    certain_infection.reserve(probs.size());

    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        p_none *= (1.0 - probs[p]);

        if (probs[p] > (1 - 1e-100))
            certain_infection.push_back(p);
        
    }

    TDbl r = static_cast<TDbl>(m->runif());
    // If there are one or more probs that go close to 1, sample
    // uniformly
    if (certain_infection.size() > 0)
        return certain_infection[std::floor(r * certain_infection.size())];

    // Step 2: Calculating the prob of none or single
    std::vector< TDbl > probs_only_p(probs.size());
    TDbl p_none_or_single = p_none;
    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        probs_only_p[p] = probs[p] * (p_none / (1.0 - probs[p]));
        p_none_or_single += probs_only_p[p];
    }

    // Step 3: Roulette
    TDbl cumsum = p_none/p_none_or_single;
    if (r < cumsum)
    {
        return -1;
    }

    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        // If it yield here, then bingo, the individual will acquire the disease
        cumsum += probs_only_p[p]/(p_none_or_single);
        if (r < cumsum)
            return static_cast<int>(p);
        
    }


    #ifdef EPI_DEBUG
    printf_epiworld("[epi-debug] roulette::cumsum = %.4f\n", cumsum);
    #endif

    return static_cast<int>(probs.size() - 1u);

}

template<typename TSeq>
inline int roulette(std::vector< double > & probs, Model<TSeq> * m)
{
    return roulette<TSeq, double>(probs, m);
}

template<typename TSeq>
inline int roulette(std::vector< float > & probs, Model<TSeq> * m)
{
    return roulette<TSeq, float>(probs, m);
}


template<typename TSeq>
inline int roulette(
    epiworld_fast_uint nelements,
    Model<TSeq> * m
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

    // Step 1: Computing the prob on none 
    epiworld_double p_none = 1.0;
    epiworld_fast_uint ncertain = 0u;
    // std::vector< int > certain_infection;
    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
    {
        p_none *= (1.0 - m->array_double_tmp[p]);

        if (m->array_double_tmp[p] > (1 - 1e-100))
            m->array_double_tmp[nelements + ncertain++] = p;
            // certain_infection.push_back(p);
        
    }

    epiworld_double r = m->runif();
    // If there are one or more probs that go close to 1, sample
    // uniformly
    if (ncertain > 0u)
        return m->array_double_tmp[nelements + std::floor(ncertain * r)]; //    certain_infection[std::floor(r * certain_infection.size())];

    // Step 2: Calculating the prob of none or single
    // std::vector< epiworld_double > probs_only_p;
    epiworld_double p_none_or_single = p_none;
    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
    {
        m->array_double_tmp[nelements + p] = 
            m->array_double_tmp[p] * (p_none / (1.0 - m->array_double_tmp[p]));
        p_none_or_single += m->array_double_tmp[nelements + p];
    }

    // Step 3: Roulette
    epiworld_double cumsum = p_none/p_none_or_single;
    if (r < cumsum)
        return -1;

    for (epiworld_fast_uint p = 0u; p < nelements; ++p)
    {
        // If it yield here, then bingo, the individual will acquire the disease
        cumsum += m->array_double_tmp[nelements + p]/(p_none_or_single);
        if (r < cumsum)
            return static_cast<int>(p);
        
    }

    return static_cast<int>(nelements - 1u);

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

#endif
