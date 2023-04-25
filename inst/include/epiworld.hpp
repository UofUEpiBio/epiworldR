#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>
#include <random>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <chrono>
#include <climits>
#include <cstdint>
#include <algorithm>
#include <regex>

#ifndef EPIWORLD_HPP
#define EPIWORLD_HPP

namespace epiworld {

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/config.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_CONFIG_HPP
#define EPIWORLD_CONFIG_HPP

#ifndef printf_epiworld
    #define printf_epiworld fflush(stdout);printf
#endif

#ifndef EPIWORLD_MAXNEIGHBORS
    #define EPIWORLD_MAXNEIGHBORS 100000
#endif

#ifdef _OPENMP
    #include <omp.h>
// #else
//     #define omp_get_thread_num() 0
//     #define omp_set_num_threads() 1
#endif

#ifndef epiworld_double
    #define epiworld_double float
#endif

#ifndef epiworld_fast_int
    #define epiworld_fast_int long long int
#endif

#ifndef epiworld_fast_uint
    #define epiworld_fast_uint unsigned long long int
#endif

#define EPI_DEFAULT_TSEQ int

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Model;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Agent;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class PersonTools;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Virus;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Tool;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Entity;

template<typename TSeq>
using VirusPtr = std::shared_ptr< Virus< TSeq > >;

template<typename TSeq>
using ToolPtr = std::shared_ptr< Tool< TSeq > >;

template<typename TSeq>
using ToolFun = std::function<epiworld_double(Tool<TSeq>&,Agent<TSeq>*,VirusPtr<TSeq>,Model<TSeq>*)>;

template<typename TSeq>
using MixerFun = std::function<epiworld_double(Agent<TSeq>*,VirusPtr<TSeq>,Model<TSeq>*)>;

template<typename TSeq>
using MutFun = std::function<bool(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq>
using PostRecoveryFun = std::function<void(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq>
using VirusFun = std::function<epiworld_double(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq>
using UpdateFun = std::function<void(Agent<TSeq>*,Model<TSeq>*)>;

template<typename TSeq>
using GlobalFun = std::function<void(Model<TSeq>*)>;

template<typename TSeq>
struct Action;

template<typename TSeq>
using ActionFun = std::function<void(Action<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute viruses at initialization
 */
template<typename TSeq>
using VirusToAgentFun = std::function<void(Virus<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute tools at initialization
 */
template<typename TSeq>
using ToolToAgentFun = std::function<void(Tool<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute entities at initialization
 */
template<typename TSeq>
using EntityToAgentFun = std::function<void(Entity<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Action data for update an agent
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
struct Action {
    Agent<TSeq> * agent;
    VirusPtr<TSeq> virus;
    ToolPtr<TSeq> tool;
    Entity<TSeq> * entity;
    epiworld_fast_int new_state;
    epiworld_fast_int queue;
    ActionFun<TSeq> call;
    int idx_agent;
    int idx_object;
public:
/**
     * @brief Construct a new Action object
     * 
     * All the parameters are rather optional.
     * 
     * @param agent_ Agent over who the action will happen
     * @param virus_ Virus to add
     * @param tool_ Tool to add
     * @param virus_idx Index of virus to be removed (if needed)
     * @param tool_idx Index of tool to be removed (if needed)
     * @param new_state_ Next state
     * @param queue_ Efect on the queue
     * @param call_ The action call (if needed)
     * @param idx_agent_ Location of agent in object.
     * @param idx_object_ Location of object in agent.
     */
    Action(
        Agent<TSeq> * agent_,
        VirusPtr<TSeq> virus_,
        ToolPtr<TSeq> tool_,
        Entity<TSeq> * entity_,
        epiworld_fast_int new_state_,
        epiworld_fast_int queue_,
        ActionFun<TSeq> call_,
        int idx_agent_,
        int idx_object_
    ) : agent(agent_), virus(virus_), tool(tool_), entity(entity_),
        new_state(new_state_),
        queue(queue_), call(call_), idx_agent(idx_agent_), idx_object(idx_object_) {
            return;
        };
};

/**
 * @name Constants in epiworld 
 * 
 * @details The following are the default values some probabilities and
 * rates take when no value has been specified in the model.
 */
///@{
#ifndef DEFAULT_TOOL_CONTAGION_REDUCTION
    #define DEFAULT_TOOL_CONTAGION_REDUCTION    0.0
#endif

#ifndef DEFAULT_TOOL_TRANSMISSION_REDUCTION
    #define DEFAULT_TOOL_TRANSMISSION_REDUCTION 0.0
#endif

#ifndef DEFAULT_TOOL_RECOVERY_ENHANCER
    #define DEFAULT_TOOL_RECOVERY_ENHANCER      0.0
#endif

#ifndef DEFAULT_TOOL_DEATH_REDUCTION
    #define DEFAULT_TOOL_DEATH_REDUCTION        0.0
#endif

#ifndef EPI_DEFAULT_VIRUS_PROB_INFECTION
    #define EPI_DEFAULT_VIRUS_PROB_INFECTION    1.0
#endif

#ifndef EPI_DEFAULT_VIRUS_PROB_RECOVERY
    #define EPI_DEFAULT_VIRUS_PROB_RECOVERY     0.1428
#endif

#ifndef EPI_DEFAULT_VIRUS_PROB_DEATH
    #define EPI_DEFAULT_VIRUS_PROB_DEATH        0.0
#endif
///@}

#ifdef EPI_DEBUG
    #define EPI_DEBUG_PRINTF printf_epiworld

    #define EPI_DEBUG_ERROR(etype, msg) \
        (etype)("[[epi-debug]] (error) " + std::string(msg));

    #define EPI_DEBUG_NOTIFY_ACTIVE() \
        EPI_DEBUG_PRINTF("DEBUGGING ON (compiled with EPI_DEBUG defined)%s\n", "");

    #define EPI_DEBUG_ALL_NON_NEGATIVE(vect) \
        for (auto & v : vect) \
            if (static_cast<double>(v) < 0.0) \
                throw EPI_DEBUG_ERROR(std::logic_error, "A negative value not allowed.");

    #define EPI_DEBUG_SUM_DBL(vect, num) \
        double _epi_debug_sum = 0.0; \
        for (auto & v : vect) \
        {   \
            _epi_debug_sum += static_cast<double>(v);\
            if (_epi_debug_sum > static_cast<double>(num)) \
                throw EPI_DEBUG_ERROR(std::logic_error, "The sum of elements not reached."); \
        }

    #define EPI_DEBUG_SUM_INT(vect, num) \
        int _epi_debug_sum = 0; \
        for (auto & v : vect) \
        {   \
            _epi_debug_sum += static_cast<int>(v);\
            if (_epi_debug_sum > static_cast<int>(num)) \
                throw EPI_DEBUG_ERROR(std::logic_error, "The sum of elements not reached."); \
        }

    #define EPI_DEBUG_VECTOR_MATCH_INT(a, b, c) \
        if (a.size() != b.size())  {\
            EPI_DEBUG_PRINTF("In '%s'", std::string(c).c_str()); \
            EPI_DEBUG_PRINTF("Size of vector a: %lu\n", (a).size());\
            EPI_DEBUG_PRINTF("Size of vector b: %lu\n", (b).size());\
            throw EPI_DEBUG_ERROR(std::length_error, "The vectors do not match size."); \
        }\
        for (int _i = 0; _i < static_cast<int>(a.size()); ++_i) \
            if (a[_i] != b[_i]) {\
                EPI_DEBUG_PRINTF("In '%s'", std::string(c).c_str()); \
                EPI_DEBUG_PRINTF("Iterating the last 5 values%s:\n", ""); \
                for (int _j = std::max(0, static_cast<int>(_i) - 4); _j <= _i; ++_j) \
                { \
                    EPI_DEBUG_PRINTF( \
                        "a[%i]: %i; b[%i]: %i\n", \
                        _j, \
                        static_cast<int>(a[_j]), \
                        _j, static_cast<int>(b[_j])); \
                } \
                throw EPI_DEBUG_ERROR(std::logic_error, "The vectors do not match."); \
            }

    #define EPI_DEBUG_FAIL_AT_TRUE(a,b) \
        if (a) \
        {\
            throw EPI_DEBUG_ERROR(std::logic_error, b); \
        } 
#else
    #define EPI_DEBUG_PRINTF(fmt, ...)
    #define EPI_DEBUG_ERROR(fmt, ...)
    #define EPI_DEBUG_NOTIFY_ACTIVE()
    #define EPI_DEBUG_ALL_NON_NEGATIVE(vect)
    #define EPI_DEBUG_SUM_DBL(vect, num)
    #define EPI_DEBUG_SUM_INT(vect, num)
    #define EPI_DEBUG_VECTOR_MATCH_INT(a, b, c)
    #define EPI_DEBUG_FAIL_AT_TRUE(a, b) \
        if (a) \
            return false;
#endif

#ifdef EPI_DEBUG_NO_THREAD_ID
    #define EPI_GET_THREAD_ID() 0
#else
    #define EPI_GET_THREAD_ID() omp_get_thread_num()
#endif

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/config.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/epiworld-macros.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MACROS_HPP
#define EPIWORLD_MACROS_HPP



/**
 * @brief Helper macro to define a new tool
 * 
 */
#define EPI_NEW_TOOL(fname,tseq) inline epiworld_double \
(fname)(\
    epiworld::Tool< tseq > & t, \
    epiworld::Agent< tseq > * p, \
    std::shared_ptr<epiworld::Virus< tseq >> v, \
    epiworld::Model< tseq > * m\
    )

/**
 * @brief Create a Tool within a function
 * 
 */
#define EPI_NEW_TOOL_LAMBDA(funname,tseq) \
    epiworld::ToolFun<tseq> funname = \
    [](epiworld::Tool<tseq> & t, \
    epiworld::Agent<tseq> * p, \
    std::shared_ptr<epiworld::Virus<tseq>> v, \
    epiworld::Model<tseq> * m) -> epiworld_double

/**
 * @brief Helper macro for accessing model parameters
 * 
 */
#define EPI_PARAMS(i) m->operator()(i)

/**
 * @brief Helper macro for defining Mutation Functions
 * 
 */
#define EPI_NEW_MUTFUN(funname,tseq) inline bool \
    (funname)(\
    epiworld::Agent<tseq> * p, \
    epiworld::Virus<tseq> & v, \
    epiworld::Model<tseq> * m )

#define EPI_NEW_MUTFUN_LAMBDA(funname,tseq) \
    epiworld::MutFun<tseq> funname = \
    [](epiworld::Agent<tseq> * p, \
    epiworld::Virus<tseq> & v, \
    epiworld::Model<tseq> * m) -> void

#define EPI_NEW_POSTRECOVERYFUN(funname,tseq) inline void \
    (funname)( \
    epiworld::Agent<tseq> * p, \
    epiworld::Virus<tseq> & v, \
    epiworld::Model<tseq> * m\
    )

#define EPI_NEW_POSTRECOVERYFUN_LAMBDA(funname,tseq) \
    epiworld::PostRecoveryFun<tseq> funname = \
    [](epiworld::Agent<tseq> * p, \
    epiworld::Virus<tseq> & v , \
    epiworld::Model<tseq> * m) -> void

#define EPI_NEW_VIRUSFUN(funname,tseq) inline epiworld_double \
    (funname)( \
    epiworld::Agent<tseq> * p, \
    epiworld::Virus<tseq> & v, \
    epiworld::Model<tseq> * m\
    )

#define EPI_NEW_VIRUSFUN_LAMBDA(funname,TSeq) \
    epiworld::VirusFun<TSeq> funname = \
    [](epiworld::Agent<TSeq> * p, \
    epiworld::Virus<TSeq> & v, \
    epiworld::Model<TSeq> * m) -> epiworld_double

#define EPI_RUNIF() m->runif()

#define EPIWORLD_RUN(a) \
    if (a.get_verbose()) \
    { \
        printf_epiworld("Running the model...\n");\
    } \
    for (epiworld_fast_uint niter = 0; niter < a.get_ndays(); ++niter)

#define EPI_TOKENPASTE(a,b) a ## b
#define MPAR(num) *(m->EPI_TOKENPASTE(p,num))

#define EPI_NEW_UPDATEFUN(funname,tseq) inline void \
    (funname)(epiworld::Agent<tseq> * p, epiworld::Model<tseq> * m)

#define EPI_NEW_UPDATEFUN_LAMBDA(funname,tseq) \
    epiworld::UpdateFun<tseq> funname = \
    [](epiworld::Agent<tseq> * p, epiworld::Model<tseq> * m) -> void

#define EPI_NEW_GLOBALFUN(funname,tseq) inline void \
    (funname)(epiworld::Model<tseq>* m)

#define EPI_NEW_GLOBALFUN_LAMBDA(funname,tseq) \
    epiworld::GlobalFun<tseq> funname = \
    [](epiworld::Model<tseq>* m) -> void

class QueueValues {
public:
    static const int NoOne    = 0;
    static const int OnlySelf = 1;
    static const int Everyone = 2;
};


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/epiworld-macros.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/misc.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


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
template<typename TSeq = int>
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
template<typename TSeq>
inline int roulette(
    const std::vector< epiworld_double > & probs,
    Model<TSeq> * m
    )
{

    // Step 1: Computing the prob on none 
    epiworld_double p_none = 1.0;
    std::vector< int > certain_infection;
    certain_infection.reserve(probs.size());

    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        p_none *= (1.0 - probs[p]);

        if (probs[p] > (1 - 1e-100))
            certain_infection.push_back(p);
        
    }

    epiworld_double r = m->runif();
    // If there are one or more probs that go close to 1, sample
    // uniformly
    if (certain_infection.size() > 0)
        return certain_infection[std::floor(r * certain_infection.size())];

    // Step 2: Calculating the prob of none or single
    std::vector< epiworld_double > probs_only_p(probs.size());
    epiworld_double p_none_or_single = p_none;
    for (epiworld_fast_uint p = 0u; p < probs.size(); ++p)
    {
        probs_only_p[p] = probs[p] * (p_none / (1.0 - probs[p]));
        p_none_or_single += probs_only_p[p];
    }

    // Step 3: Roulette
    epiworld_double cumsum = p_none/p_none_or_single;
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
inline int roulette(
    epiworld_fast_uint nelements,
    Model<TSeq> * m
    )
{

    #ifdef EPI_DEBUG
    if (nelements > m->array_double_tmp.size())
        throw std::logic_error("Trying to sample from more data than there is in roulette!");
    #endif

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

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/misc.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/progress.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_PROGRESS_HPP
#define EPIWORLD_PROGRESS_HPP

#ifndef EPIWORLD_PROGRESS_BAR_WIDTH
#define EPIWORLD_PROGRESS_BAR_WIDTH 80
#endif

/**
 * @brief A simple progress bar
  */
class Progress {
private:
    int    width;     ///< Total width size (number of bars)
    int    n;         ///< Total number of iterations
    epiworld_double step_size; ///< Size of the step
    int last_loc;     ///< Last location of the bar
    int cur_loc;      ///< Last location of the bar
    int i;            ///< Current iteration step
    
public:
    Progress() {};
    Progress(int n_, int width_);
    ~Progress() {};

    void start();
    void next();
    void end();

};

inline Progress::Progress(int n_, int width_) {


    width     = std::max(7, width_ - 7);
    n         = n_;
    step_size = static_cast<epiworld_double>(width)/static_cast<epiworld_double>(n);
    last_loc  = 0;
    i         = 0;

}

inline void Progress::start()
{

    #ifndef EPI_DEBUG
    for (int j = 0; j < (width); ++j)
    {
        printf_epiworld("_");
    }
    printf_epiworld("\n");
    #endif
}

inline void Progress::next() {

    if (i == 0)
        start();

    cur_loc = std::floor((++i) * step_size);


    #ifndef EPI_DEBUG
    for (int j = 0; j < (cur_loc - last_loc); ++j)
    {
        printf_epiworld("|");
    }
    #endif

    if (i >= n)
        end();

    last_loc = cur_loc;

}

inline void Progress::end() {

    #ifndef EPI_DEBUG
    printf_epiworld(" done.\n");
    #endif

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/progress.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



    // #include "math/summary-stats.hpp"

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/math/lfmcmc.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_LFMCMC_HPP
#define EPIWORLD_LFMCMC_HPP

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//math/lfmcmc/lfmcmc-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_LFMCMC_BONES_HPP
#define EPIWORLD_LFMCMC_BONES_HPP

#ifndef epiworld_double
    #define epiworld_double float
#endif


template<typename TData>
class LFMCMC;

template<typename TData>
using LFMCMCSimFun = std::function<TData(const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCSummaryFun = std::function<void(std::vector< epiworld_double >&,const TData&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCProposalFun = std::function<void(std::vector< epiworld_double >&,const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCKernelFun = std::function<epiworld_double(const std::vector< epiworld_double >&,const std::vector< epiworld_double >&,epiworld_double,LFMCMC<TData>*)>;

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
);

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
    epiworld_double lb = std::numeric_limits<epiworld_double>::min(),
    epiworld_double ub = std::numeric_limits<epiworld_double>::max()
);

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
);

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
);

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
);

/**
 * @brief Likelihood-Free Markov Chain Monte Carlo
 * 
 * @tparam TData Type of data that is generated
 */
template<typename TData>
class LFMCMC {
private:

    // Random number sampling
    std::mt19937 * engine = nullptr;
    
    std::shared_ptr< std::uniform_real_distribution<> > runifd =
        std::make_shared< std::uniform_real_distribution<> >(0.0, 1.0);

    std::shared_ptr< std::normal_distribution<> > rnormd =
        std::make_shared< std::normal_distribution<> >(0.0);

    std::shared_ptr< std::gamma_distribution<> > rgammad = 
        std::make_shared< std::gamma_distribution<> >();

    // Process data
    TData * observed_data;
    
    // Information about the size of the problem
    size_t n_samples;
    size_t n_statistics;
    size_t n_parameters;

    epiworld_double epsilon;

    std::vector< epiworld_double > params_now;
    std::vector< epiworld_double > params_prev;
    std::vector< epiworld_double > params_init;

    std::vector< epiworld_double > observed_stats; ///< Observed statistics

    std::vector< epiworld_double > sampled_params;     ///< Sampled Parameters
    std::vector< epiworld_double > sampled_stats;      ///< Sampled statistics
    std::vector< epiworld_double > sampled_stats_prob; ///< Sampled statistics
    std::vector< bool >            sampled_accepted;   ///< Indicator of accepted statistics

    std::vector< epiworld_double > accepted_params;      ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_stats;       ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_params_prob; ///< Posterior probability

    std::vector< epiworld_double > drawn_prob;     ///< Drawn probabilities (runif())
    std::vector< TData > * sampled_data = nullptr;

    // Functions
    LFMCMCSimFun<TData> simulation_fun;
    LFMCMCSummaryFun<TData> summary_fun;
    LFMCMCProposalFun<TData> proposal_fun = proposal_fun_normal<TData>;
    LFMCMCKernelFun<TData> kernel_fun     = kernel_fun_uniform<TData>;

    // Misc
    std::vector< std::string > names_parameters;
    std::vector< std::string > names_statistics;

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed = 
        std::chrono::duration<epiworld_double,std::micro>::zero();

    inline void get_elapsed(
        std::string unit,
        epiworld_double * last_elapsed,
        std::string * unit_abbr,
        bool print
    );

    void chrono_start();
    void chrono_end();
    
public:

    void run(
        std::vector< epiworld_double > param_init,
        size_t n_samples_,
        epiworld_double epsilon_
        );

    LFMCMC() {};
    LFMCMC(TData & observed_data_) : observed_data(&observed_data_) {};
    ~LFMCMC() {};

    void set_observed_data(TData & observed_data_) {observed_data = &observed_data_;};
    void set_proposal_fun(LFMCMCProposalFun<TData> fun);
    void set_simulation_fun(LFMCMCSimFun<TData> fun);
    void set_summary_fun(LFMCMCSummaryFun<TData> fun);
    void set_kernel_fun(LFMCMCKernelFun<TData> fun);
    
    /**
     * @name Random number generation
     * 
     * @param eng 
     */
    ///@{
    void set_rand_engine(std::mt19937 & eng);
    std::mt19937 & get_rand_endgine();
    void seed(epiworld_fast_uint s);
    void set_rand_gamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double runif();
    epiworld_double rnorm();
    epiworld_double rgamma();
    epiworld_double runif(epiworld_double lb, epiworld_double ub);
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    ///@}

    // Accessing parameters of the function
    size_t get_n_samples() const {return n_samples;};
    size_t get_n_statistics() const {return n_statistics;};
    size_t get_n_parameters() const {return n_parameters;};
    epiworld_double get_epsilon() const {return epsilon;};

    const std::vector< epiworld_double > & get_params_now() {return params_now;};
    const std::vector< epiworld_double > & get_params_prev() {return params_prev;};
    const std::vector< epiworld_double > & get_params_init() {return params_init;};
    const std::vector< epiworld_double > & get_statistics_obs() {return observed_stats;};
    const std::vector< epiworld_double > & get_statistics_hist() {return sampled_stats;};
    const std::vector< bool >            & get_statistics_accepted() {return sampled_accepted;};
    const std::vector< epiworld_double > & get_posterior_lf_prob() {return accepted_params_prob;};
    const std::vector< epiworld_double > & get_drawn_prob() {return drawn_prob;};
    std::vector< TData > * get_sampled_data() {return sampled_data;};

    void set_par_names(std::vector< std::string > names);
    void set_stats_names(std::vector< std::string > names);

    std::vector< epiworld_double > get_params_mean();
    std::vector< epiworld_double > get_stats_mean();

    void print() ;

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//math/lfmcmc/lfmcmc-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//math/lfmcmc/lfmcmc-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_LFMCMC_MEAT_HPP
#define EPIWORLD_LFMCMC_MEAT_HPP

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//math//lfmcmc/lfmcmc-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_LFMCMC_BONES_HPP
#define EPIWORLD_LFMCMC_BONES_HPP

#ifndef epiworld_double
    #define epiworld_double float
#endif


template<typename TData>
class LFMCMC;

template<typename TData>
using LFMCMCSimFun = std::function<TData(const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCSummaryFun = std::function<void(std::vector< epiworld_double >&,const TData&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCProposalFun = std::function<void(std::vector< epiworld_double >&,const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCKernelFun = std::function<epiworld_double(const std::vector< epiworld_double >&,const std::vector< epiworld_double >&,epiworld_double,LFMCMC<TData>*)>;

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
);

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
    epiworld_double lb = std::numeric_limits<epiworld_double>::min(),
    epiworld_double ub = std::numeric_limits<epiworld_double>::max()
);

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
);

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
);

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
);

/**
 * @brief Likelihood-Free Markov Chain Monte Carlo
 * 
 * @tparam TData Type of data that is generated
 */
template<typename TData>
class LFMCMC {
private:

    // Random number sampling
    std::mt19937 * engine = nullptr;
    
    std::shared_ptr< std::uniform_real_distribution<> > runifd =
        std::make_shared< std::uniform_real_distribution<> >(0.0, 1.0);

    std::shared_ptr< std::normal_distribution<> > rnormd =
        std::make_shared< std::normal_distribution<> >(0.0);

    std::shared_ptr< std::gamma_distribution<> > rgammad = 
        std::make_shared< std::gamma_distribution<> >();

    // Process data
    TData * observed_data;
    
    // Information about the size of the problem
    size_t n_samples;
    size_t n_statistics;
    size_t n_parameters;

    epiworld_double epsilon;

    std::vector< epiworld_double > params_now;
    std::vector< epiworld_double > params_prev;
    std::vector< epiworld_double > params_init;

    std::vector< epiworld_double > observed_stats; ///< Observed statistics

    std::vector< epiworld_double > sampled_params;     ///< Sampled Parameters
    std::vector< epiworld_double > sampled_stats;      ///< Sampled statistics
    std::vector< epiworld_double > sampled_stats_prob; ///< Sampled statistics
    std::vector< bool >            sampled_accepted;   ///< Indicator of accepted statistics

    std::vector< epiworld_double > accepted_params;      ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_stats;       ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_params_prob; ///< Posterior probability

    std::vector< epiworld_double > drawn_prob;     ///< Drawn probabilities (runif())
    std::vector< TData > * sampled_data = nullptr;

    // Functions
    LFMCMCSimFun<TData> simulation_fun;
    LFMCMCSummaryFun<TData> summary_fun;
    LFMCMCProposalFun<TData> proposal_fun = proposal_fun_normal<TData>;
    LFMCMCKernelFun<TData> kernel_fun     = kernel_fun_uniform<TData>;

    // Misc
    std::vector< std::string > names_parameters;
    std::vector< std::string > names_statistics;

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed = 
        std::chrono::duration<epiworld_double,std::micro>::zero();

    inline void get_elapsed(
        std::string unit,
        epiworld_double * last_elapsed,
        std::string * unit_abbr,
        bool print
    );

    void chrono_start();
    void chrono_end();
    
public:

    void run(
        std::vector< epiworld_double > param_init,
        size_t n_samples_,
        epiworld_double epsilon_
        );

    LFMCMC() {};
    LFMCMC(TData & observed_data_) : observed_data(&observed_data_) {};
    ~LFMCMC() {};

    void set_observed_data(TData & observed_data_) {observed_data = &observed_data_;};
    void set_proposal_fun(LFMCMCProposalFun<TData> fun);
    void set_simulation_fun(LFMCMCSimFun<TData> fun);
    void set_summary_fun(LFMCMCSummaryFun<TData> fun);
    void set_kernel_fun(LFMCMCKernelFun<TData> fun);
    
    /**
     * @name Random number generation
     * 
     * @param eng 
     */
    ///@{
    void set_rand_engine(std::mt19937 & eng);
    std::mt19937 & get_rand_endgine();
    void seed(epiworld_fast_uint s);
    void set_rand_gamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double runif();
    epiworld_double rnorm();
    epiworld_double rgamma();
    epiworld_double runif(epiworld_double lb, epiworld_double ub);
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    ///@}

    // Accessing parameters of the function
    size_t get_n_samples() const {return n_samples;};
    size_t get_n_statistics() const {return n_statistics;};
    size_t get_n_parameters() const {return n_parameters;};
    epiworld_double get_epsilon() const {return epsilon;};

    const std::vector< epiworld_double > & get_params_now() {return params_now;};
    const std::vector< epiworld_double > & get_params_prev() {return params_prev;};
    const std::vector< epiworld_double > & get_params_init() {return params_init;};
    const std::vector< epiworld_double > & get_statistics_obs() {return observed_stats;};
    const std::vector< epiworld_double > & get_statistics_hist() {return sampled_stats;};
    const std::vector< bool >            & get_statistics_accepted() {return sampled_accepted;};
    const std::vector< epiworld_double > & get_posterior_lf_prob() {return accepted_params_prob;};
    const std::vector< epiworld_double > & get_drawn_prob() {return drawn_prob;};
    std::vector< TData > * get_sampled_data() {return sampled_data;};

    void set_par_names(std::vector< std::string > names);
    void set_stats_names(std::vector< std::string > names);

    std::vector< epiworld_double > get_params_mean();
    std::vector< epiworld_double > get_stats_mean();

    void print() ;

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//math//lfmcmc/lfmcmc-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



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
    summary_fun(observed_stats, *observed_data, this);
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
inline void LFMCMC<TData>::set_rand_engine(std::mt19937 & eng)
{
    engine = &eng;
}

template<typename TData>
inline void LFMCMC<TData>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::make_shared<std::gamma_distribution<>>(alpha,beta);
}

template<typename TData>
inline std::mt19937 & LFMCMC<TData>::get_rand_endgine()
{
    return *engine;
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

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//math//lfmcmc/lfmcmc-meat-print.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef LFMCMC_MEAT_PRINT_HPP
#define LFMCMC_MEAT_PRINT_HPP

template<typename TData>
inline void LFMCMC<TData>::print()
{
    std::vector< epiworld_double > summ_params(n_parameters * 3, 0.0);
    std::vector< epiworld_double > summ_stats(n_statistics * 3, 0.0);

    for (size_t k = 0u; k < n_parameters; ++k)
    {

        // Retrieving the relevant parameter
        std::vector< epiworld_double > par_i(n_samples);
        for (size_t i = 0u; i < n_samples; ++i)
        {
            par_i[i] = accepted_params[i * n_parameters + k];
            summ_params[k * 3] += par_i[i]/n_samples;
        }

        // Computing the 95% Credible interval
        std::sort(par_i.begin(), par_i.end());

        summ_params[k * 3 + 1u] = 
            par_i[std::floor(.025 * static_cast<epiworld_double>(n_samples))];
        summ_params[k * 3 + 2u] = 
            par_i[std::floor(.975 * static_cast<epiworld_double>(n_samples))];

    }

    for (size_t k = 0u; k < n_statistics; ++k)
    {

        // Retrieving the relevant parameter
        std::vector< epiworld_double > stat_k(n_samples);
        for (size_t i = 0u; i < n_samples; ++i)
        {
            stat_k[i] = accepted_stats[i * n_statistics + k];
            summ_stats[k * 3] += stat_k[i]/n_samples;
        }

        // Computing the 95% Credible interval
        std::sort(stat_k.begin(), stat_k.end());

        summ_stats[k * 3 + 1u] = 
            stat_k[std::floor(.025 * static_cast<epiworld_double>(n_samples))];
        summ_stats[k * 3 + 2u] = 
            stat_k[std::floor(.975 * static_cast<epiworld_double>(n_samples))];

    }

    printf_epiworld("___________________________________________\n\n");
    printf_epiworld("LIKELIHOOD-FREE MARKOV CHAIN MONTE CARLO\n\n");

    printf_epiworld("N Samples : %ld\n", n_samples);

    std::string abbr;
    epiworld_double elapsed;
    get_elapsed("auto", &elapsed, &abbr, false);
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
        
        int tmp_nchar = std::floor(std::log10(std::abs(n)));
        if (nchar_par_num < tmp_nchar)
            nchar_par_num = tmp_nchar;
    }
    nchar_par_num += 5; // 1 for neg padd, 2 for decimals, 1 the decimal point, and one b/c log(<10) < 1.
    std::string charlen = std::to_string(nchar_par_num);

    if (names_parameters.size() != 0u)
    {
        int nchar_par = 0;
        for (auto & n : names_parameters)
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

        for (size_t k = 0u; k < n_parameters; ++k)
        {
            printf_epiworld(
                fmt_params.c_str(),
                names_parameters[k].c_str(),
                summ_params[k * 3],
                summ_params[k * 3 + 1u],
                summ_params[k * 3 + 2u],
                params_init[k]
                );
        }

        
    } else {

        fmt_params = std::string("  [%-2ld]: % ") + charlen + 
            std::string(".2f [% ") + charlen +
            std::string(".2f, % ") + charlen + 
            std::string(".2f] (initial : % ") + charlen +
            std::string(".2f)\n");

        for (size_t k = 0u; k < n_parameters; ++k)
        {
            
            printf_epiworld(
                fmt_params.c_str(),
                k,
                summ_params[k * 3],
                summ_params[k * 3 + 1u],
                summ_params[k * 3 + 2u],
                params_init[k]
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
        int tmp_nchar = std::floor(std::log10(std::abs(s)));
        if (nchar < tmp_nchar)
            nchar = tmp_nchar;
    }

    nchar += 5; // See above

    std::string nchar_char = std::to_string(nchar);

    // Figuring out format
    std::string fmt_stats;
    if (names_statistics.size() != 0u)
    {
        int nchar_stats = 0;
        for (auto & n : names_statistics)
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

        for (size_t k = 0u; k < n_statistics; ++k)
        {
            printf_epiworld(
                fmt_stats.c_str(),
                names_statistics[k].c_str(),
                summ_stats[k * 3],
                summ_stats[k * 3 + 1u],
                summ_stats[k * 3 + 2u],
                observed_stats[k]
                );
        }

        
    } else {

        fmt_stats = std::string("  [%-2ld] : % ") + 
            nchar_char +
            std::string(".2f [% ") + nchar_char +
            std::string(".2f, % ") + nchar_char +
            std::string(".2f] (Observed: % ") + nchar_char +
            std::string(".2f)\n");

        for (size_t k = 0u; k < n_statistics; ++k)
        {
            printf_epiworld(
                fmt_stats.c_str(),
                k,
                summ_stats[k * 3],
                summ_stats[k * 3 + 1u],
                summ_stats[k * 3 + 2u],
                observed_stats[k]
                );
        }

    }

    printf_epiworld("___________________________________________\n\n");
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//math//lfmcmc/lfmcmc-meat-print.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



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
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//math/lfmcmc/lfmcmc-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/math/lfmcmc.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/userdata-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_USERDATA_BONES_HPP
#define EPIWORLD_USERDATA_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class DataBase;

/**
 * @brief Personalized data by the user
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class UserData
{
    friend class Model<TSeq>;
    friend class DataBase<TSeq>;

private:
    Model<TSeq> * model;

    std::vector< std::string > data_names;
    std::vector< int > data_dates;
    std::vector< epiworld_double > data_data;

    epiworld_fast_uint k = 0u;
    epiworld_fast_uint n = 0u;

    int last_day = -1;

public:

    UserData() = delete;
    UserData(Model<TSeq> & m) : model(&m) {};
    UserData(Model<TSeq> * m) : model(m) {};

    /**
     * @brief Construct a new User Data object
     * 
     * @param names A vector of names. The length of the vector sets
     * the number of columns to record.
     */
    UserData(std::vector< std::string > names);

    /**
     * @name Append data 
     * 
     * @param x A vector of length `ncol()` (if vector), otherwise a `epiworld_double`.
     * @param j Index of the data point, from 0 to `ncol() - 1`.
     */
    ///@{
    void add(std::vector<epiworld_double> x);
    void add(
        epiworld_fast_uint j,
        epiworld_double x
        );
    ///@}

    /**
     * @name Access data 
     * 
     * @param i Row (0 through ndays - 1.)
     * @param j Column (0 through `ncols()`).
     * @return epiworld_double& 
     */
    ///@{
    epiworld_double & operator()(
        epiworld_fast_uint i,
        epiworld_fast_uint j
        );

    epiworld_double & operator()(
        epiworld_fast_uint i,
        std::string name
        );
    ///@}

    std::vector< std::string > & get_names();

    std::vector< int > & get_dates();

    std::vector< epiworld_double > & get_data();

    void get_all(
        std::vector< std::string > * names    = nullptr,
        std::vector< int > * date             = nullptr,
        std::vector< epiworld_double > * data = nullptr
    );

    epiworld_fast_uint nrow() const;
    epiworld_fast_uint ncol() const;

    void write(std::string fn);
    void print() const;

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/userdata-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/userdata-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_USERDATA_MEAT_HPP
#define EPIWORLD_USERDATA_MEAT_HPP

template<typename TSeq>
class UserData;

template<typename TSeq>
inline UserData<TSeq>::UserData(std::vector< std::string > names)
{

    k = names.size();
    data_names = names;

}

template<typename TSeq>
inline void UserData<TSeq>::add(std::vector<epiworld_double> x)
{

    if (x.size() != k)
        throw std::out_of_range(
            "The size of -x-, " + std::to_string(x.size()) + ", does not match " +
            "the number of elements registered (" + std::to_string(k));

    for (auto & i : x)
        data_data.push_back(i);

    data_dates.push_back(model->today());

    n++;
    last_day = model->today();

}

template<typename TSeq>
inline void UserData<TSeq>::add(epiworld_fast_uint j, epiworld_double x)
{

    // Starting with a new day?
    if (static_cast<int>(model->today()) != last_day)
    {

        std::vector< epiworld_double > tmp(k, 0.0);

        tmp[j] = x;

        add(tmp);

    }
    else
    {

        this->operator()(n - 1, j) = x;

    }

}

template<typename TSeq>
inline std::vector< std::string > & UserData<TSeq>::get_names() 
{
    return data_names;
}

template<typename TSeq>
inline std::vector< int > & UserData<TSeq>::get_dates() 
{
    return data_dates;
}

template<typename TSeq>
inline std::vector< epiworld_double > & UserData<TSeq>::get_data() 
{
    return data_data;
}

template<typename TSeq>
inline void UserData<TSeq>::get_all(
    std::vector< std::string > * names,
    std::vector< int > * dates,
    std::vector< epiworld_double > * data
) 
{
    
    if (names != nullptr)
        names = &this->data_names;

    if (dates != nullptr)
        dates = &this->data_dates;

    if (data != nullptr)
        data = &this->data_data;

}

template<typename TSeq>
inline epiworld_double & UserData<TSeq>::operator()(
    epiworld_fast_uint i,
    epiworld_fast_uint j
)
{

    if (j >= k)
        throw std::out_of_range("j cannot be greater than k - 1.");

    if (i >= n)
        throw std::out_of_range("j cannot be greater than n - 1.");

    return data_data[k * i + j];

}

template<typename TSeq>
inline epiworld_double & UserData<TSeq>::operator()(
    epiworld_fast_uint i,
    std::string name
)
{
    int loc = -1;
    for (epiworld_fast_uint l = 0u; l < k; ++l)
    {

        if (name == data_names[l])
        {

            loc = l;
            break;

        }

    }

    if (loc < 0)
        throw std::range_error(
            "The variable \"" + name + "\" is not present " +
            "in the user UserData database."
        );

    return operator()(i, static_cast<epiworld_fast_uint>(loc));

}

template<typename TSeq>
inline epiworld_fast_uint UserData<TSeq>::nrow() const
{
    return n;
}

template<typename TSeq>
inline epiworld_fast_uint UserData<TSeq>::ncol() const
{
    return k;
}

template<typename TSeq>
inline void UserData<TSeq>::write(std::string fn)
{
    std::ofstream file_ud(fn, std::ios_base::out);

    // File header
    file_ud << "\"date\"";
    for (auto & cn : data_names)
        file_ud << " \"" + cn + "\"";
    file_ud << "\n";
    
    epiworld_fast_uint ndata = 0u;
    for (epiworld_fast_uint i = 0u; i < n; ++i)
    {
        file_ud << data_dates[i];

        for (epiworld_fast_uint j = 0u; j < k; ++j)
            file_ud << " " << data_data[ndata++];

        file_ud << "\n";
    }

    return;
}

template<typename TSeq>
inline void UserData<TSeq>::print() const
{
    // File header
    printf_epiworld("Total records: %llu\n", n);
    printf_epiworld("date");

    for (auto & cn : data_names)
    {

        printf_epiworld(" %s", cn.c_str());

    }

    printf_epiworld("\n");
    
    epiworld_fast_uint ndata = 0u;
    
    for (epiworld_fast_uint i = 0u; i < n; ++i)
    {

        printf_epiworld("%i", data_dates[i]);

        for (epiworld_fast_uint j = 0u; j < k; ++j)
        {

            printf_epiworld(" %.2f", data_data[ndata++]);

        }

        printf_epiworld("\n");

    }

    return;
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/userdata-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/seq_processing.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_SEQ_PROCESSING_HPP 
#define EPIWORLD_SEQ_PROCESSING_HPP

/**
 * @brief Hasher function to turn the sequence into an integer vector
 * 
 * @tparam TSeq 
 * @param x 
 * @return std::vector<int> 
 */
template<typename TSeq>
inline std::vector<int> default_seq_hasher(const TSeq & x);

template<>
inline std::vector<int> default_seq_hasher<std::vector<int>>(const std::vector<int> & x) {
    return x;
}

template<>
inline std::vector<int> default_seq_hasher<std::vector<bool>>(const std::vector<bool> & x) {
    std::vector<int> ans(x.size());
    size_t j = 0;
    for (const auto & i : x)
        ans[j++] = i? 1 : 0;
    return ans;
}

template<>
inline std::vector<int> default_seq_hasher<int>(const int & x) {
    return {x};
}

template<>
inline std::vector<int> default_seq_hasher<bool>(const bool & x) {
    return {x ? 1 : 0};
}

/**
 * @brief Default way to write sequences
 * 
 * @tparam TSeq 
 * @param seq 
 * @return std::string 
 */
template<typename TSeq = int>
inline std::string default_seq_writer(const TSeq & seq);

template<>
inline std::string default_seq_writer<std::vector<int>>(
    const std::vector<int> & seq
) {

    std::string out = "";
    for (const auto & s : seq)
        out = out + std::to_string(s);

    return out;

}

template<>
inline std::string default_seq_writer<std::vector<bool>>(
    const std::vector<bool> & seq
) {

    std::string out = "";
    for (const auto & s : seq)
        out = out + (s ? "1" : "0");

    return out;

}

template<>
inline std::string default_seq_writer<bool>(
    const bool & seq
) {

    return seq ? "1" : "0";

}

template<>
inline std::string default_seq_writer<int>(
    const int & seq
) {

    return std::to_string(seq);

}



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/seq_processing.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/database-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_DATABASE_BONES_HPP
#define EPIWORLD_DATABASE_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Virus;

template<typename TSeq>
class UserData;

template<typename TSeq>
inline void default_add_virus(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_tool(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_virus(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_tool(Action<TSeq> & a, Model<TSeq> * m);

/**
 * @brief Statistical data about the process
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class DataBase {
    friend class Model<TSeq>;
    friend void default_add_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_add_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:
    Model<TSeq> * model;

    // Variants information 
    MapVec_type<int,int> variant_id; ///< The squence is the key
    std::vector< std::string > variant_name;
    std::vector< TSeq> variant_sequence;
    std::vector< int > variant_origin_date;
    std::vector< int > variant_parent_id;

    MapVec_type<int,int> tool_id; ///< The squence is the key
    std::vector< std::string > tool_name;
    std::vector< TSeq> tool_sequence;
    std::vector< int > tool_origin_date;

    std::function<std::vector<int>(const TSeq&)> seq_hasher = default_seq_hasher<TSeq>;
    std::function<std::string(const TSeq &)> seq_writer = default_seq_writer<TSeq>;

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_variant;

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_tool;

    // {Susceptible, Infected, etc.}
    std::vector< int > today_total;

    // Totals
    int today_total_nvariants_active = 0;
    
    int sampling_freq = 1;

    // Variants history
    std::vector< int > hist_variant_date;
    std::vector< int > hist_variant_id;
    std::vector< epiworld_fast_uint > hist_variant_state;
    std::vector< int > hist_variant_counts;

    // Tools history
    std::vector< int > hist_tool_date;
    std::vector< int > hist_tool_id;
    std::vector< epiworld_fast_uint > hist_tool_state;
    std::vector< int > hist_tool_counts;

    // Overall hist
    std::vector< int > hist_total_date;
    std::vector< int > hist_total_nvariants_active;
    std::vector< epiworld_fast_uint > hist_total_state;
    std::vector< int > hist_total_counts;
    std::vector< int > hist_transition_matrix;

    // Transmission network
    std::vector< int > transmission_date;                 ///< Date of the transmission event
    std::vector< int > transmission_source;               ///< Id of the source
    std::vector< int > transmission_target;               ///< Id of the target
    std::vector< int > transmission_variant;              ///< Id of the variant
    std::vector< int > transmission_source_exposure_date; ///< Date when the source acquired the variant

    std::vector< int > transition_matrix;

    UserData<TSeq> user_data;

    void update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo = false
    );

    void update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
    );

    void update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
    );

    void record_transition(epiworld_fast_uint from, epiworld_fast_uint to, bool undo);


public:

    #ifdef EPI_DEBUG
    int n_transmissions_potential = 0;
    int n_transmissions_today     = 0;
    #endif

    DataBase() = delete;
    DataBase(Model<TSeq> & m) : model(&m), user_data(m) {};
    DataBase(const DataBase<TSeq> & db);
    // DataBase<TSeq> & operator=(const DataBase<TSeq> & m);

    /**
     * @brief Registering a new variant
     * 
     * @param v Pointer to the new variant.
     * Since variants are originated in the agent, the numbers simply move around.
     * From the parent variant to the new variant. And the total number of infected
     * does not change.
     */
    void record_variant(Virus<TSeq> & v); 
    void record_tool(Tool<TSeq> & t); 
    void set_seq_hasher(std::function<std::vector<int>(TSeq)> fun);
    void reset();
    Model<TSeq> * get_model();
    void record();

    const std::vector< TSeq > & get_sequence() const;
    const std::vector< int > & get_nexposed() const;
    size_t size() const;

    /**
     * @name Get recorded information from the model
     * 
     * @param what std::string, The state, e.g., 0, 1, 2, ...
     * @return In `get_today_total`, the current counts of `what`.
     * @return In `get_today_variant`, the current counts of `what` for
     * each variant.
     * @return In `get_hist_total`, the time series of `what`
     * @return In `get_hist_variant`, the time series of what for each variant.
     * @return In `get_hist_total_date` and `get_hist_variant_date` the
     * corresponding dates
     */
    ///@{
    int get_today_total(std::string what) const;
    int get_today_total(epiworld_fast_uint what) const;
    void get_today_total(
        std::vector< std::string > * state = nullptr,
        std::vector< int > * counts = nullptr
    ) const;

    void get_today_variant(
        std::vector< std::string > & state,
        std::vector< int > & id,
        std::vector< int > & counts
    ) const;

    void get_hist_total(
        std::vector< int > * date,
        std::vector< std::string > * state,
        std::vector< int > * counts
    ) const;

    void get_hist_variant(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & state,
        std::vector< int > & counts
    ) const;

    void get_hist_tool(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & state,
        std::vector< int > & counts
    ) const;

    void get_hist_transition_matrix(
        std::vector< std::string > & state_from,
        std::vector< std::string > & state_to,
        std::vector< int > & date,
        std::vector< int > & counts,
        bool skip_zeros
    ) const;
    ///@}

    void write_data(
        std::string fn_variant_info,
        std::string fn_variant_hist,
        std::string fn_tool_info,
        std::string fn_tool_hist,
        std::string fn_total_hist,
        std::string fn_transmission,
        std::string fn_transition,
        std::string fn_reproductive_number,
        std::string fn_generation_time
        ) const;
    
    void record_transmission(int i, int j, int variant, int i_expo_date);

    size_t get_n_variants() const;
    size_t get_n_tools() const;
    
    void set_user_data(std::vector< std::string > names);
    void add_user_data(std::vector< epiworld_double > x);
    void add_user_data(epiworld_fast_uint j, epiworld_double x);
    UserData<TSeq> & get_user_data();


    /**
     * @brief Computes the reproductive number of each case
     * 
     * @details By definition, whereas it computes R0 (basic reproductive number)
     * or Rt/R (the effective reproductive number) will depend on whether the
     * virus is allowed to circulate naïvely or not, respectively.
     * 
     * @param fn File where to write out the reproductive number.
     */
    ///@{
    MapVec_type<int,int> reproductive_number() const;

    void reproductive_number(
        std::string fn
        ) const;
    ///@}

    /**
     * @brief Calculates the transition probabilities
     * 
     * @return std::vector< epiworld_double > 
     */
    std::vector< epiworld_double > transition_probability(
        bool print = true
    ) const;

    bool operator==(const DataBase<TSeq> & other) const;
    bool operator!=(const DataBase<TSeq> & other) const {return !operator==(other);};

    /**
     * Calculates the generating time
     * @param agent_id,virus_id,time,gentime vectors where to save the values agent_id
    */
   ///@{
    void generation_time(
        std::vector< int > & agent_id,
        std::vector< int > & virus_id,
        std::vector< int > & time,
        std::vector< int > & gentime
    ) const;

    void generation_time(
        std::string fn
    ) const;
    ///@}

};


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/database-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/database-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_DATABASE_MEAT_HPP
#define EPIWORLD_DATABASE_MEAT_HPP

template<typename TSeq>
inline void DataBase<TSeq>::reset()
{

    // Initializing the counts
    today_total.resize(model->nstatus);
    std::fill(today_total.begin(), today_total.end(), 0);
    for (auto & p : model->get_agents())
        ++today_total[p.get_state()];
    
    transition_matrix.resize(model->nstatus * model->nstatus);
    std::fill(transition_matrix.begin(), transition_matrix.end(), 0);
    for (size_t s = 0u; s < model->nstatus; ++s)
        transition_matrix[s + s * model->nstatus] = today_total[s];

    hist_variant_date.clear();
    hist_variant_id.clear();
    hist_variant_state.clear();
    hist_variant_counts.clear();

    hist_tool_date.clear();
    hist_tool_id.clear();
    hist_tool_state.clear();
    hist_tool_counts.clear();    

    today_variant.resize(get_n_variants());
    std::fill(today_variant.begin(), today_variant.begin(), std::vector<int>(model->nstatus, 0));

    today_tool.resize(get_n_tools());
    std::fill(today_tool.begin(), today_tool.begin(), std::vector<int>(model->nstatus, 0));

    hist_total_date.clear();
    hist_total_state.clear();
    hist_total_nvariants_active.clear();
    hist_total_counts.clear();
    hist_transition_matrix.clear();

    transmission_date.clear();
    transmission_variant.clear();
    transmission_source.clear();
    transmission_target.clear();
    transmission_source_exposure_date.clear();

    

    return;

}

template<typename TSeq>
inline DataBase<TSeq>::DataBase(const DataBase<TSeq> & db) :
    variant_id(db.variant_id),
    variant_name(db.variant_name),
    variant_sequence(db.variant_sequence),
    variant_origin_date(db.variant_origin_date),
    variant_parent_id(db.variant_parent_id),
    tool_id(db.tool_id),
    tool_name(db.tool_name),
    tool_sequence(db.tool_sequence),
    tool_origin_date(db.tool_origin_date),
    seq_hasher(db.seq_hasher),
    seq_writer(db.seq_writer),
    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    today_variant(db.today_variant),
    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    today_tool(db.today_tool),
    // {Susceptible, Infected, etc.}
    today_total(db.today_total),
    // Totals
    today_total_nvariants_active(db.today_total_nvariants_active),
    sampling_freq(db.sampling_freq),
    // Variants history
    hist_variant_date(db.hist_variant_date),
    hist_variant_id(db.hist_variant_id),
    hist_variant_state(db.hist_variant_state),
    hist_variant_counts(db.hist_variant_counts),
    // Tools history
    hist_tool_date(db.hist_tool_date),
    hist_tool_id(db.hist_tool_id),
    hist_tool_state(db.hist_tool_state),
    hist_tool_counts(db.hist_tool_counts),
    // Overall hist
    hist_total_date(db.hist_total_date),
    hist_total_nvariants_active(db.hist_total_nvariants_active),
    hist_total_state(db.hist_total_state),
    hist_total_counts(db.hist_total_counts),
    hist_transition_matrix(db.hist_transition_matrix),
    // Transmission network
    transmission_date(db.transmission_date),
    transmission_source(db.transmission_source),
    transmission_target(db.transmission_target),
    transmission_variant(db.transmission_variant),
    transmission_source_exposure_date(db.transmission_source_exposure_date),
    transition_matrix(db.transition_matrix),
    user_data(nullptr)
{}

// DataBase<TSeq> & DataBase<TSeq>::operator=(const DataBase<TSeq> & m)
// {

// }

template<typename TSeq>
inline Model<TSeq> * DataBase<TSeq>::get_model() {
    return model;
}

template<typename TSeq>
inline const std::vector< TSeq > & DataBase<TSeq>::get_sequence() const {
    return variant_sequence;
}

template<typename TSeq>
inline void DataBase<TSeq>::record() 
{

    ////////////////////////////////////////////////////////////////////////////
    // DEBUGGING BLOCK
    ////////////////////////////////////////////////////////////////////////////
    EPI_DEBUG_SUM_INT(today_total, model->size())
    EPI_DEBUG_ALL_NON_NEGATIVE(today_total)

    #ifdef EPI_DEBUG
    // Checking whether the sums correspond
    std::vector< int > _today_total_cp(today_total.size(), 0);
    for (auto & p : model->population)
        _today_total_cp[p.get_state()]++;
    
    EPI_DEBUG_VECTOR_MATCH_INT(
        _today_total_cp, today_total,
        "Sums of __today_total_cp in database-meat.hpp"
        )

    if (model->today() == 0)
    {
        if (hist_total_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_date should be of length 0.")
        if (hist_total_nvariants_active.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_nvariants_active should be of length 0.")
        if (hist_total_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_state should be of length 0.")
        if (hist_total_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_counts should be of length 0.")
        if (hist_variant_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_variant_date should be of length 0.")
        if (hist_variant_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_variant_id should be of length 0.")
        if (hist_variant_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_variant_state should be of length 0.")
        if (hist_variant_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_variant_counts should be of length 0.")
        if (hist_tool_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_date should be of length 0.")
        if (hist_tool_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_id should be of length 0.")
        if (hist_tool_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_state should be of length 0.")
        if (hist_tool_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_counts should be of length 0.")
    }
    #endif
    ////////////////////////////////////////////////////////////////////////////

    // Only store every now and then
    if ((model->today() % sampling_freq) == 0)
    {

        // Recording variant's history
        for (auto & p : variant_id)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstatus; ++s)
            {

                hist_variant_date.push_back(model->today());
                hist_variant_id.push_back(p.second);
                hist_variant_state.push_back(s);
                hist_variant_counts.push_back(today_variant[p.second][s]);

            }

        }

        // Recording tool's history
        for (auto & p : tool_id)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstatus; ++s)
            {

                hist_tool_date.push_back(model->today());
                hist_tool_id.push_back(p.second);
                hist_tool_state.push_back(s);
                hist_tool_counts.push_back(today_tool[p.second][s]);

            }

        }

        // Recording the overall history
        for (epiworld_fast_uint s = 0u; s < model->nstatus; ++s)
        {
            hist_total_date.push_back(model->today());
            hist_total_nvariants_active.push_back(today_total_nvariants_active);
            hist_total_state.push_back(s);
            hist_total_counts.push_back(today_total[s]);
        }

        for (auto cell : transition_matrix)
            hist_transition_matrix.push_back(cell);

        // Now the diagonal must reflect the state
        for (size_t s_i = 0u; s_i < model->nstatus; ++s_i)
        {
            for (size_t s_j = 0u; s_j < model->nstatus; ++s_j)
            {
                if ((s_i != s_j) && (transition_matrix[s_i + s_j * model->nstatus] > 0))
                {
                    transition_matrix[s_j + s_j * model->nstatus] +=
                        transition_matrix[s_i + s_j * model->nstatus];

                    transition_matrix[s_i + s_j * model->nstatus] = 0;
                }
         
            }

            #ifdef EPI_DEBUG
            if (transition_matrix[s_i + s_i * model->nstatus] != 
                today_total[s_i])
                throw std::logic_error(
                    "The diagonal of the updated transition Matrix should match the daily totals"
                    );
            #endif
        }


    }

}

template<typename TSeq>
inline void DataBase<TSeq>::record_variant(Virus<TSeq> & v)
{

    // If no sequence, then need to add one. This is regardless of the case
    if (v.get_sequence() == nullptr)
        v.set_sequence(default_sequence<TSeq>(
            static_cast<int>(variant_name.size())
            ));

    // Negative id -> virus hasn't been recorded
    if (v.get_id() < 0)
    {


        // Generating the hash
        std::vector< int > hash = seq_hasher(*v.get_sequence());

        epiworld_fast_uint new_id = variant_id.size();
        variant_id[hash] = new_id;
        variant_name.push_back(v.get_name());
        variant_sequence.push_back(*v.get_sequence());
        variant_origin_date.push_back(model->today());
        
        variant_parent_id.push_back(v.get_id()); // Must be -99
        
        today_variant.push_back({});
        today_variant[new_id].resize(model->nstatus, 0);
       
        // Updating the variant
        v.set_id(new_id);
        v.set_date(model->today());

        today_total_nvariants_active++;

    } else { // In this case, the virus is already on record, need to make sure
             // The new sequence is new.

        // Updating registry
        std::vector< int > hash = seq_hasher(*v.get_sequence());
        epiworld_fast_uint old_id = v.get_id();
        epiworld_fast_uint new_id;

        // If the sequence is new, then it means that the
        if (variant_id.find(hash) == variant_id.end())
        {

            new_id = variant_id.size();
            variant_id[hash] = new_id;
            variant_name.push_back(v.get_name());
            variant_sequence.push_back(*v.get_sequence());
            variant_origin_date.push_back(model->today());
            
            variant_parent_id.push_back(old_id);
            
            today_variant.push_back({});
            today_variant[new_id].resize(model->nstatus, 0);
        
            // Updating the variant
            v.set_id(new_id);
            v.set_date(model->today());

            today_total_nvariants_active++;

        } else {

            // Finding the id
            new_id = variant_id[hash];

            // Reflecting the change
            v.set_id(new_id);
            v.set_date(variant_origin_date[new_id]);

        }

        // Moving statistics (only if we are affecting an individual)
        if (v.get_agent() != nullptr)
        {
            // Correcting math
            epiworld_fast_uint tmp_state = v.get_agent()->get_state();
            today_variant[old_id][tmp_state]--;
            today_variant[new_id][tmp_state]++;

        }

    }
    
    return;

} 

template<typename TSeq>
inline void DataBase<TSeq>::record_tool(Tool<TSeq> & t)
{

    if (t.get_sequence() == nullptr)
        t.set_sequence(default_sequence<TSeq>(
            static_cast<int>(tool_name.size())
        ));

    if (t.get_id() < 0) 
    {

        std::vector< int > hash = seq_hasher(*t.get_sequence());
        epiworld_fast_uint new_id = tool_id.size();
        tool_id[hash] = new_id;
        tool_name.push_back(t.get_name());
        tool_sequence.push_back(*t.get_sequence());
        tool_origin_date.push_back(model->today());
                
        today_tool.push_back({});
        today_tool[new_id].resize(model->nstatus, 0);

        // Updating the tool
        t.set_id(new_id);
        t.set_date(model->today());

    } else {

        // Updating registry
        std::vector< int > hash = seq_hasher(*t.get_sequence());
        epiworld_fast_uint old_id = t.get_id();
        epiworld_fast_uint new_id;
        
        if (tool_id.find(hash) == tool_id.end())
        {

            new_id = tool_id.size();
            tool_id[hash] = new_id;
            tool_name.push_back(t.get_name());
            tool_sequence.push_back(*t.get_sequence());
            tool_origin_date.push_back(model->today());
                    
            today_tool.push_back({});
            today_tool[new_id].resize(model->nstatus, 0);

            // Updating the tool
            t.set_id(new_id);
            t.set_date(model->today());

        } else {

            // Finding the id
            new_id = tool_id[hash];

            // Reflecting the change
            t.set_id(new_id);
            t.set_date(tool_origin_date[new_id]);

        }

        // Moving statistics (only if we are affecting an individual)
        if (t.get_agent() != nullptr)
        {
            // Correcting math
            epiworld_fast_uint tmp_state = t.get_agent()->get_state();
            today_tool[old_id][tmp_state]--;
            today_tool[new_id][tmp_state]++;

        }

    }

    
    
    return;
} 

template<typename TSeq>
inline size_t DataBase<TSeq>::size() const
{
    return variant_id.size();
}

template<typename TSeq>
inline void DataBase<TSeq>::update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    if (undo)
    {

        today_total[prev_state]++;
        today_total[new_state]--;
        
    } else {

        today_total[prev_state]--;
        today_total[new_state]++;

    }

    record_transition(prev_state, new_state, undo);
    
    return;
}

template<typename TSeq>
inline void DataBase<TSeq>::update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
) {

    today_variant[virus_id][prev_state]--;
    today_variant[virus_id][new_state]++;

    return;
    
}

template<typename TSeq>
inline void DataBase<TSeq>::update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
) {


    today_tool[tool_id][prev_state]--;    
    today_tool[tool_id][new_state]++;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::record_transition(
    epiworld_fast_uint from,
    epiworld_fast_uint to,
    bool undo
) {

    if (undo)
    {   

        transition_matrix[to * model->nstatus + from]--;
        transition_matrix[from * model->nstatus + from]++;

    } else {

        transition_matrix[to * model->nstatus + from]++;
        transition_matrix[from * model->nstatus + from]--;

    }

    #ifdef EPI_DEBUG
    if (transition_matrix[from * model->nstatus + from] < 0)
        throw std::logic_error("An entry in transition matrix is negative.");
    #endif

}

template<typename TSeq>
inline int DataBase<TSeq>::get_today_total(
    std::string what
) const
{

    for (auto i = 0u; i < model->states_labels.size(); ++i)
    {
        if (model->states_labels[i] == what)
            return today_total[i];
    }

    throw std::range_error("The value '" + what + "' is not in the model.");

}

template<typename TSeq>
inline void DataBase<TSeq>::get_today_total(
    std::vector< std::string > * state,
    std::vector< int > * counts
) const
{
    if (state != nullptr)
        (*state) = model->states_labels;

    if (counts != nullptr)
        *counts = today_total;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_today_variant(
    std::vector< std::string > & state,
    std::vector< int > & id,
    std::vector< int > & counts
    ) const
{
      
    state.resize(today_variant.size(), "");
    id.resize(today_variant.size(), 0);
    counts.resize(today_variant.size(),0);

    int n = 0u;
    for (epiworld_fast_uint v = 0u; v < today_variant.size(); ++v)
        for (epiworld_fast_uint s = 0u; s < model->states_labels.size(); ++s)
        {
            state[n]   = model->states_labels[s];
            id[n]       = static_cast<int>(v);
            counts[n++] = today_variant[v][s];

        }

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_total(
    std::vector< int > * date,
    std::vector< std::string > * state,
    std::vector< int > * counts
) const
{

    if (date != nullptr)
        *date = hist_total_date;

    if (state != nullptr)
    {
        state->resize(hist_total_state.size(), "");
        for (epiworld_fast_uint i = 0u; i < hist_total_state.size(); ++i)
            state->operator[](i) = model->states_labels[hist_total_state[i]];
    }

    if (counts != nullptr)
        *counts = hist_total_counts;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_variant(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & state,
    std::vector< int > & counts
) const {

    date = hist_variant_date;
    std::vector< std::string > labels;
    labels = model->states_labels;
    
    id = hist_variant_id;
    state.resize(hist_variant_state.size(), "");
    for (epiworld_fast_uint i = 0u; i < hist_variant_state.size(); ++i)
        state[i] = labels[hist_variant_state[i]];

    counts = hist_variant_counts;

    return;

}


template<typename TSeq>
inline void DataBase<TSeq>::get_hist_tool(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & state,
    std::vector< int > & counts
) const {

    date = hist_tool_date;
    std::vector< std::string > labels;
    labels = model->states_labels;
    
    id = hist_tool_id;
    state.resize(hist_tool_state.size(), "");
    for (size_t i = 0u; i < hist_tool_state.size(); ++i)
        state[i] = labels[hist_tool_state[i]];

    counts = hist_tool_counts;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_transition_matrix(
    std::vector< std::string > & state_from,
    std::vector< std::string > & state_to,
    std::vector< int > & date,
    std::vector< int > & counts,
    bool skip_zeros
) const
{

    size_t n = this->hist_transition_matrix.size();
    
    // Reserving space
    state_from.reserve(n);
    state_to.reserve(n);
    date.reserve(n);
    counts.reserve(n);

    size_t n_status = model->nstatus;
    size_t n_steps  = model->get_ndays();

    for (size_t step = 0u; step <= n_steps; ++step) // The final step counts
    {
        for (size_t j = 0u; j < n_status; ++j) // Column major storage
        {
            for (size_t i = 0u; i < n_status; ++i)
            {
                // Retrieving the value of the day
                int v = hist_transition_matrix[
                    step * n_status * n_status + // Day of the data
                    j * n_status +               // Column (to)
                    i                            // Row (from)
                    ];

                // If we are skipping the zeros and it is zero, then don't save
                if (skip_zeros && v == 0)
                    continue;
                                
                state_from.push_back(model->states_labels[i]);
                state_to.push_back(model->states_labels[j]);
                date.push_back(hist_total_date[step * n_status]);
                counts.push_back(v);

            }

        }
    }

    return;


}

template<typename TSeq>
inline void DataBase<TSeq>::write_data(
    std::string fn_variant_info,
    std::string fn_variant_hist,
    std::string fn_tool_info,
    std::string fn_tool_hist,
    std::string fn_total_hist,
    std::string fn_transmission,
    std::string fn_transition,
    std::string fn_reproductive_number,
    std::string fn_generation_time
) const
{

    if (fn_variant_info != "")
    {
        std::ofstream file_variant_info(fn_variant_info, std::ios_base::out);

        file_variant_info <<
        #ifdef _OPENMP
            "thread" << "id " << "variant_name " << "variant_sequence " << "date_recorded " << "parent\n";
        #else
            "id " << "variant_name " << "variant_sequence " << "date_recorded " << "parent\n";
        #endif

        for (const auto & v : variant_id)
        {
            int id = v.second;
            file_variant_info <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                id << " \"" <<
                variant_name[id] << "\" " <<
                seq_writer(variant_sequence[id]) << " " <<
                variant_origin_date[id] << " " <<
                variant_parent_id[id] << "\n";
        }

    }

    if (fn_variant_hist != "")
    {
        std::ofstream file_variant(fn_variant_hist, std::ios_base::out);
        
        file_variant <<
            #ifdef _OPENMP
            "thread "<< "date " << "id " << "state " << "n\n";
            #else
            "date " << "id " << "state " << "n\n";
            #endif

        for (epiworld_fast_uint i = 0; i < hist_variant_id.size(); ++i)
            file_variant <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_variant_date[i] << " " <<
                hist_variant_id[i] << " " <<
                model->states_labels[hist_variant_state[i]] << " " <<
                hist_variant_counts[i] << "\n";
    }

    if (fn_tool_info != "")
    {
        std::ofstream file_tool_info(fn_tool_info, std::ios_base::out);

        file_tool_info <<
            #ifdef _OPENMP
            "thread " << 
            #endif
            "id " << "tool_name " << "tool_sequence " << "date_recorded\n";

        for (const auto & t : tool_id)
        {
            int id = t.second;
            file_tool_info <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                id << " \"" <<
                tool_name[id] << "\" " <<
                seq_writer(tool_sequence[id]) << " " <<
                tool_origin_date[id] << "\n";
        }

    }

    if (fn_tool_hist != "")
    {
        std::ofstream file_tool_hist(fn_tool_hist, std::ios_base::out);
        
        file_tool_hist <<
            #ifdef _OPENMP
            "thread " << 
            #endif
            "date " << "id " << "state " << "n\n";

        for (epiworld_fast_uint i = 0; i < hist_tool_id.size(); ++i)
            file_tool_hist <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_tool_date[i] << " " <<
                hist_tool_id[i] << " " <<
                model->states_labels[hist_tool_state[i]] << " " <<
                hist_tool_counts[i] << "\n";
    }

    if (fn_total_hist != "")
    {
        std::ofstream file_total(fn_total_hist, std::ios_base::out);

        file_total <<
            #ifdef _OPENMP
            "thread " << 
            #endif
            "date " << "nvariants " << "state " << "counts\n";

        for (epiworld_fast_uint i = 0; i < hist_total_date.size(); ++i)
            file_total <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_total_date[i] << " " <<
                hist_total_nvariants_active[i] << " \"" <<
                model->states_labels[hist_total_state[i]] << "\" " << 
                hist_total_counts[i] << "\n";
    }

    if (fn_transmission != "")
    {
        std::ofstream file_transmission(fn_transmission, std::ios_base::out);
        file_transmission <<
            #ifdef _OPENMP
            "thread " << 
            #endif
            "date " << "variant " << "source_exposure_date " << "source " << "target\n";

        for (epiworld_fast_uint i = 0; i < transmission_target.size(); ++i)
            file_transmission <<
                #ifdef _OPENMP
                EPI_GET_THREAD_ID() << " " <<
                #endif
                transmission_date[i] << " " <<
                transmission_variant[i] << " " <<
                transmission_source_exposure_date[i] << " " <<
                transmission_source[i] << " " <<
                transmission_target[i] << "\n";
                
    }

    if (fn_transition != "")
    {
        std::ofstream file_transition(fn_transition, std::ios_base::out);
        file_transition <<
            #ifdef _OPENMP
            "thread " << 
            #endif
            "date " << "from " << "to " << "counts\n";

        int ns = model->nstatus;

        for (int i = 0; i <= model->today(); ++i)
        {

            for (int from = 0u; from < ns; ++from)
                for (int to = 0u; to < ns; ++to)
                    file_transition <<
                        #ifdef _OPENMP
                        EPI_GET_THREAD_ID() << " " <<
                        #endif
                        i << " " <<
                        model->states_labels[from] << " " <<
                        model->states_labels[to] << " " <<
                        hist_transition_matrix[i * (ns * ns) + to * ns + from] << "\n";
                
        }
                
    }

    if (fn_reproductive_number != "")
        reproductive_number(fn_reproductive_number);

    if (fn_generation_time != "")
        generation_time(fn_generation_time);

}

template<typename TSeq>
inline void DataBase<TSeq>::record_transmission(
    int i,
    int j,
    int variant,
    int i_expo_date
) {

    transmission_date.push_back(model->today());
    transmission_source.push_back(i);
    transmission_target.push_back(j);
    transmission_variant.push_back(variant);
    transmission_source_exposure_date.push_back(i_expo_date);

}

template<typename TSeq>
inline size_t DataBase<TSeq>::get_n_variants() const
{
    return variant_id.size();
}

template<typename TSeq>
inline size_t DataBase<TSeq>::get_n_tools() const
{
    return tool_id.size();
}


template<typename TSeq>
inline void DataBase<TSeq>::set_user_data(
    std::vector< std::string > names
)
{
    user_data = UserData<TSeq>(names);
    user_data.model = model;
}

template<typename TSeq>
inline void DataBase<TSeq>::add_user_data(
    std::vector< epiworld_double > x
)
{

    user_data.add(x);

}

template<typename TSeq>
inline void DataBase<TSeq>::add_user_data(
    epiworld_fast_uint k,
    epiworld_double x
)
{

    user_data.add(k, x);

}

template<typename TSeq>
inline UserData<TSeq> & DataBase<TSeq>::get_user_data()
{
    return user_data;
}

template<typename TSeq>
inline MapVec_type<int,int> DataBase<TSeq>::reproductive_number()
const {

    // Checking size
    MapVec_type<int,int> map;

    // Number of digits of maxid
    for (size_t i = 0u; i < transmission_date.size(); ++i)
    {
        // Fabricating id
        std::vector< int > h = {
            transmission_variant[i],
            transmission_source[i],
            transmission_source_exposure_date[i]
        };

        // Adding to counter
        if (map.find(h) == map.end())
            map[h] = 1;
        else
            map[h]++;

        // The target is added
        std::vector< int > h_target = {
            transmission_variant[i],
            transmission_target[i],
            transmission_date[i]
        };

        map[h_target] = 0;
        
    }

    return map;

}

template<typename TSeq>
inline void DataBase<TSeq>::reproductive_number(
    std::string fn
) const {


    auto map = reproductive_number();

    std::ofstream fn_file(fn, std::ios_base::out);

    fn_file << 
        #ifdef _OPENMP
        "thread " <<
        #endif
        "variant source source_exposure_date rt\n";

    for (auto & m : map)
        fn_file <<
            #ifdef _OPENMP
            EPI_GET_THREAD_ID() << " " <<
            #endif
            m.first[0u] << " " <<
            m.first[1u] << " " <<
            m.first[2u] << " " <<
            m.second << "\n";

    return;

}

template<typename TSeq>
inline std::vector< epiworld_double > DataBase<TSeq>::transition_probability(
    bool print
) const {

    auto states_labels = model->get_state();
    size_t n_state = states_labels.size();
    size_t n_days   = model->get_ndays();
    std::vector< epiworld_double > res(n_state * n_state, 0.0);
    std::vector< epiworld_double > days_to_include(n_state, 0.0);

    for (size_t t = 1; t < n_days; ++t)
    {

        for (size_t s_i = 0; s_i < n_state; ++s_i)
        {
            epiworld_double daily_total = hist_total_counts[(t - 1) * n_state + s_i];

            if (daily_total == 0)
                continue;

            days_to_include[s_i] += 1.0; 

            for (size_t s_j = 0u; s_j < n_state; ++s_j)
            {
                #ifdef EPI_DEBUG
                epiworld_double entry = hist_transition_matrix[
                    s_i + s_j * n_state +
                    t * (n_state * n_state)
                    ];

                if (entry > daily_total)
                    throw std::logic_error(
                        "The entry in hist_transition_matrix cannot have more elememnts than the total"
                        );

                res[s_i + s_j * n_state] += (entry / daily_total);
                #else
                    res[s_i + s_j * n_state] += (
                        hist_transition_matrix[
                            s_i + s_j * n_state +
                            t * (n_state * n_state)
                        ] / daily_total
                    );
                #endif
            }

        }

    }

    for (size_t s_i = 0; s_i < n_state; ++s_i)
    {
        for (size_t s_j = 0; s_j < n_state; ++s_j)
            res[s_i + s_j * n_state] /= days_to_include[s_i];
    }

    if (print)
    {   

        size_t nchar = 0u;
        for (auto & l : states_labels)
            if (l.length() > nchar)
                nchar = l.length();

        std::string fmt = " - %-" + std::to_string(nchar) + "s";
        
        printf_epiworld("\nTransition Probabilities:\n");
        for (size_t s_i = 0u; s_i < n_state; ++s_i)
        {
            printf_epiworld(fmt.c_str(), states_labels[s_i].c_str());
            for (size_t s_j = 0u; s_j < n_state; ++s_j)
            {
                if (std::isnan(res[s_i + s_j * n_state]))
                {
                    printf_epiworld("     -");
                } else {
                    printf_epiworld(" % 4.2f", res[s_i + s_j * n_state]);
                }
            }
            printf_epiworld("\n");
        }

        printf_epiworld("\n");

    }

    return res;


} 

#define VECT_MATCH(a, b, c) \
    EPI_DEBUG_FAIL_AT_TRUE(a.size() != b.size(), c) \
    for (size_t __i = 0u; __i < a.size(); ++__i) \
    {\
        EPI_DEBUG_FAIL_AT_TRUE(a[__i] != b[__i], c) \
    }

template<>
inline bool DataBase<std::vector<int>>::operator==(const DataBase<std::vector<int>> & other) const
{
    VECT_MATCH(
        variant_name, other.variant_name,
        "DataBase:: variant_name don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        variant_sequence.size() != other.variant_sequence.size(),
        "DataBase:: variant_sequence don't match."
        )

    for (size_t i = 0u; i < variant_sequence.size(); ++i)
    {
        VECT_MATCH(
            variant_sequence[i], other.variant_sequence[i],
            "DataBase:: variant_sequence[i] don't match"
            )
    }

    VECT_MATCH(
        variant_origin_date,
        other.variant_origin_date,
        "DataBase:: variant_origin_date[i] don't match"
    )

    VECT_MATCH(
        variant_parent_id,
        other.variant_parent_id,
        "DataBase:: variant_parent_id[i] don't match"
    )

    VECT_MATCH(
        tool_name,
        other.tool_name,
        "DataBase:: tool_name[i] don't match"
    )

    VECT_MATCH(
        tool_sequence,
        other.tool_sequence,
        "DataBase:: tool_sequence[i] don't match"
    )

    VECT_MATCH(
        tool_origin_date,
        other.tool_origin_date,
        "DataBase:: tool_origin_date[i] don't match"
    )


    EPI_DEBUG_FAIL_AT_TRUE(
        sampling_freq != other.sampling_freq,
        "DataBase:: sampling_freq don't match."
        )

    // Variants history
    VECT_MATCH(
        hist_variant_date,
        other.hist_variant_date,
        "DataBase:: hist_variant_date[i] don't match"
        )

    VECT_MATCH(
        hist_variant_id,
        other.hist_variant_id,
        "DataBase:: hist_variant_id[i] don't match"
        )

    VECT_MATCH(
        hist_variant_state,
        other.hist_variant_state,
        "DataBase:: hist_variant_state[i] don't match"
        )

    VECT_MATCH(
        hist_variant_counts,
        other.hist_variant_counts,
        "DataBase:: hist_variant_counts[i] don't match"
        )

    // Tools history
    VECT_MATCH(
        hist_tool_date,
        other.hist_tool_date,
        "DataBase:: hist_tool_date[i] don't match"
        )

    VECT_MATCH(
        hist_tool_id,
        other.hist_tool_id,
        "DataBase:: hist_tool_id[i] don't match"
        )

    VECT_MATCH(
        hist_tool_state,
        other.hist_tool_state,
        "DataBase:: hist_tool_state[i] don't match"
        )

    VECT_MATCH(
        hist_tool_counts,
        other.hist_tool_counts,
        "DataBase:: hist_tool_counts[i] don't match"
        )

    // Overall hist
    VECT_MATCH(
        hist_total_date,
        other.hist_total_date,
        "DataBase:: hist_total_date[i] don't match"
        )

    VECT_MATCH(
        hist_total_nvariants_active,
        other.hist_total_nvariants_active,
        "DataBase:: hist_total_nvariants_active[i] don't match"
        )

    VECT_MATCH(
        hist_total_state,
        other.hist_total_state,
        "DataBase:: hist_total_state[i] don't match"
        )

    VECT_MATCH(
        hist_total_counts,
        other.hist_total_counts,
        "DataBase:: hist_total_counts[i] don't match"
        )

    VECT_MATCH(
        hist_transition_matrix,
        other.hist_transition_matrix,
        "DataBase:: hist_transition_matrix[i] don't match"
        )

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    EPI_DEBUG_FAIL_AT_TRUE(
        today_variant.size() != other.today_variant.size(),
        "DataBase:: today_variant don't match."
        )
    
    for (size_t i = 0u; i < today_variant.size(); ++i)
    {
        VECT_MATCH(
            today_variant[i], other.today_variant[i],
            "DataBase:: today_variant[i] don't match"
            )
    }

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    if (today_tool.size() != other.today_tool.size())
        return false;
    
    for (size_t i = 0u; i < today_tool.size(); ++i)
    {
        VECT_MATCH(
            today_tool[i], other.today_tool[i],
            "DataBase:: today_tool[i] don't match"
            )
    }

    // {Susceptible, Infected, etc.}
    VECT_MATCH(
        today_total, other.today_total,
        "DataBase:: today_total don't match"
        )

    // Totals
    EPI_DEBUG_FAIL_AT_TRUE(
        today_total_nvariants_active != other.today_total_nvariants_active,
        "DataBase:: today_total_nvariants_active don't match."
        )

    // Transmission network
    VECT_MATCH(
        transmission_date,
        other.transmission_date,                 ///< Date of the transmission eve,
        "DataBase:: transmission_date[i] don't match"
        )

    VECT_MATCH(
        transmission_source,
        other.transmission_source,               ///< Id of the sour,
        "DataBase:: transmission_source[i] don't match"
        )

    VECT_MATCH(
        transmission_target,
        other.transmission_target,               ///< Id of the targ,
        "DataBase:: transmission_target[i] don't match"
        )

    VECT_MATCH(
        transmission_variant,
        other.transmission_variant,              ///< Id of the varia,
        "DataBase:: transmission_variant[i] don't match"
        )

    VECT_MATCH(
        transmission_source_exposure_date,
        other.transmission_source_exposure_date, ///< Date when the source acquired the varia,
        "DataBase:: transmission_source_exposure_date[i] don't match"
        )


    VECT_MATCH(
        transition_matrix,
        other.transition_matrix,
        "DataBase:: transition_matrix[i] don't match"
        )


    return true;

}

template<typename TSeq>
inline bool DataBase<TSeq>::operator==(const DataBase<TSeq> & other) const
{
    VECT_MATCH(
        variant_name,
        other.variant_name,
        "DataBase:: variant_name[i] don't match"
    )

    VECT_MATCH(
        variant_sequence,
        other.variant_sequence,
        "DataBase:: variant_sequence[i] don't match"
    )

    VECT_MATCH(
        variant_origin_date,
        other.variant_origin_date,
        "DataBase:: variant_origin_date[i] don't match"
    )

    VECT_MATCH(
        variant_parent_id,
        other.variant_parent_id,
        "DataBase:: variant_parent_id[i] don't match"
    )

    VECT_MATCH(
        tool_name,
        other.tool_name,
        "DataBase:: tool_name[i] don't match"
    )

    VECT_MATCH(
        tool_sequence,
        other.tool_sequence,
        "DataBase:: tool_sequence[i] don't match"
    )

    VECT_MATCH(
        tool_origin_date,
        other.tool_origin_date,
        "DataBase:: tool_origin_date[i] don't match"
    )

    
    EPI_DEBUG_FAIL_AT_TRUE(
        sampling_freq != other.sampling_freq,
        "DataBase:: sampling_freq don't match."
    )

    // Variants history
    VECT_MATCH(
        hist_variant_date,
        other.hist_variant_date,
        "DataBase:: hist_variant_date[i] don't match"
    )

    VECT_MATCH(
        hist_variant_id,
        other.hist_variant_id,
        "DataBase:: hist_variant_id[i] don't match"
    )

    VECT_MATCH(
        hist_variant_state,
        other.hist_variant_state,
        "DataBase:: hist_variant_state[i] don't match"
    )

    VECT_MATCH(
        hist_variant_counts,
        other.hist_variant_counts,
        "DataBase:: hist_variant_counts[i] don't match"
    )

    // Tools history
    VECT_MATCH(
        hist_tool_date,
        other.hist_tool_date,
        "DataBase:: hist_tool_date[i] don't match"
    )

    VECT_MATCH(
        hist_tool_id,
        other.hist_tool_id,
        "DataBase:: hist_tool_id[i] don't match"
    )

    VECT_MATCH(
        hist_tool_state,
        other.hist_tool_state,
        "DataBase:: hist_tool_state[i] don't match"
    )

    VECT_MATCH(
        hist_tool_counts,
        other.hist_tool_counts,
        "DataBase:: hist_tool_counts[i] don't match"
    )

    // Overall hist
    VECT_MATCH(
        hist_total_date,
        other.hist_total_date,
        "DataBase:: hist_total_date[i] don't match"
    )

    VECT_MATCH(
        hist_total_nvariants_active,
        other.hist_total_nvariants_active,
        "DataBase:: hist_total_nvariants_active[i] don't match"
    )

    VECT_MATCH(
        hist_total_state,
        other.hist_total_state,
        "DataBase:: hist_total_state[i] don't match"
    )

    VECT_MATCH(
        hist_total_counts,
        other.hist_total_counts,
        "DataBase:: hist_total_counts[i] don't match"
    )

    VECT_MATCH(
        hist_transition_matrix,
        other.hist_transition_matrix,
        "DataBase:: hist_transition_matrix[i] don't match"
    )

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    EPI_DEBUG_FAIL_AT_TRUE(
        today_variant.size() != other.today_variant.size(),
        "DataBase:: today_variant.size() don't match."
    )
    
    for (size_t i = 0u; i < today_variant.size(); ++i)
    {
        VECT_MATCH(
            today_variant[i], other.today_variant[i],
            "DataBase:: today_variant[i] don't match"
            )
    }

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    EPI_DEBUG_FAIL_AT_TRUE(
        today_tool.size() != other.today_tool.size(),
        "DataBase:: today_tool.size() don't match."
    )
    
    for (size_t i = 0u; i < today_tool.size(); ++i)
    {
        VECT_MATCH(
            today_tool[i], other.today_tool[i],
            "DataBase:: today_tool[i] don't match"
            )
    }

    // {Susceptible, Infected, etc.}
    VECT_MATCH(
        today_total, other.today_total,
        "DataBase:: today_total[i] don't match"
        )

    // Totals
    EPI_DEBUG_FAIL_AT_TRUE(
        today_total_nvariants_active != other.today_total_nvariants_active,
        "DataBase:: today_total_nvariants_active don't match."
    )

    // Transmission network
    VECT_MATCH( ///< Date of the transmission eve
        transmission_date,
        other.transmission_date,
        "DataBase:: transmission_date[i] don't match"
    )

    VECT_MATCH( ///< Id of the sour
        transmission_source,
        other.transmission_source,
        "DataBase:: transmission_source[i] don't match"
    )

    VECT_MATCH( ///< Id of the targ
        transmission_target,
        other.transmission_target,
        "DataBase:: transmission_target[i] don't match"
    )

    VECT_MATCH( ///< Id of the varia
        transmission_variant,
        other.transmission_variant,
        "DataBase:: transmission_variant[i] don't match"
    )

    VECT_MATCH( ///< Date when the source acquired the varia
        transmission_source_exposure_date,
        other.transmission_source_exposure_date,
        "DataBase:: transmission_source_exposure_date[i] don't match"
    )

    VECT_MATCH(
        transition_matrix,
        other.transition_matrix,
        "DataBase:: transition_matrix[i] don't match"
    )

    return true;

}

template<typename TSeq>
inline void DataBase<TSeq>::generation_time(
    std::vector< int > & agent_id,
    std::vector< int > & virus_id,
    std::vector< int > & time,
    std::vector< int > & gentime
) const {
    
    size_t nevents = transmission_date.size();

    agent_id.reserve(nevents);
    virus_id.reserve(nevents);
    time.reserve(nevents);
    gentime.reserve(nevents);

    // Iterating through the individuals
    for (size_t i = 0u; i < nevents; ++i)
    {
        int agent_id_i = transmission_target[i];
        agent_id.push_back(agent_id_i);
        virus_id.push_back(transmission_variant[i]);
        time.push_back(transmission_date[i]);

        bool found = false;
        for (size_t j = i; j < nevents; ++j)
        {

            if (transmission_source[j] == agent_id_i)
            {
                gentime.push_back(transmission_date[j] - time[i]);
                found = true;
                break;
            }

        }

        // If there's no transmission, we set the generation time to
        // minus 1;
        if (!found)
            gentime.push_back(-1);

    }

    agent_id.shrink_to_fit();
    virus_id.shrink_to_fit();
    time.shrink_to_fit();
    gentime.shrink_to_fit();

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::generation_time(
    std::string fn
) const
{

    std::vector< int > agent_id;
    std::vector< int > virus_id;
    std::vector< int > time;
    std::vector< int > gentime;

    generation_time(agent_id, virus_id, time, gentime);

    std::ofstream fn_file(fn, std::ios_base::out);

    fn_file << 
        #ifdef _OPENMP
        "thread " <<
        #endif
        "variant source source_exposure_date gentime\n";

    size_t n = agent_id.size();
    for (size_t i = 0u; i < n; ++i)
        fn_file <<
            #ifdef _OPENMP
            EPI_GET_THREAD_ID() << " " <<
            #endif
            virus_id[i] << " " <<
            agent_id[i] << " " <<
            time[i] << " " <<
            gentime[i] << "\n";

    return;

}

#undef VECT_MATCH

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/database-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/adjlist-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_ADJLIST_BONES_HPP
#define EPIWORLD_ADJLIST_BONES_HPP

class AdjList {
private:

    std::vector<std::map<int, int>> dat;
    bool directed;
    epiworld_fast_uint N = 0;
    epiworld_fast_uint E = 0;

public:

    AdjList() {};

    /**
     * @brief Construct a new Adj List object
     * 
     * @details 
     * Ids in the network are assume to range from `0` to `size - 1`.
     * 
     * @param source Unsigned int vector with the source
     * @param target Unsigned int vector with the target
     * @param size Number of vertices in the network.
     * @param directed Bool true if the network is directed
     */
    AdjList(
        const std::vector< int > & source,
        const std::vector< int > & target,
        int size,
        bool directed
        );

    AdjList(AdjList && a); // Move constructor
    AdjList(const AdjList & a); // Copy constructor
    AdjList& operator=(const AdjList& a);


    /**
     * @brief Read an edgelist
     * 
     * Ids in the network are assume to range from `0` to `size - 1`.
     * 
     * @param fn Path to the file
     * @param skip Number of lines to skip (e.g., 1 if there's a header)
     * @param directed `true` if the network is directed
     * @param size Number of vertices in the network.
     */
    void read_edgelist(
        std::string fn,
        int size,
        int skip = 0,
        bool directed = true
        );

    std::map<int, int> operator()(
        epiworld_fast_uint i
        ) const;
        
    void print(epiworld_fast_uint limit = 20u) const;
    size_t vcount() const; ///< Number of vertices/nodes in the network.
    size_t ecount() const; ///< Number of edges/arcs/ties in the network.
    
    std::vector<std::map<int,int>> & get_dat() {
        return dat;
    };

    bool is_directed() const; ///< `true` if the network is directed.

};


#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/adjlist-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/adjlist-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_ADJLIST_MEAT_HPP
#define EPIWORLD_ADJLIST_MEAT_HPP

inline AdjList::AdjList(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) : directed(directed) {


    dat.resize(size, std::map<int,int>({}));
    int max_id = size - 1;

    int i,j;
    for (int m = 0; m < static_cast<int>(source.size()); ++m)
    {

        i = source[m];
        j = target[m];

        if (i > max_id)
            throw std::range_error(
                "The source["+std::to_string(m)+"] = " + std::to_string(i) +
                " is above the max_id " + std::to_string(max_id)
                );

        if (j > max_id)
            throw std::range_error(
                "The target["+std::to_string(m)+"] = " + std::to_string(j) +
                " is above the max_id " + std::to_string(max_id)
                );

        // Adding nodes
        if (dat[i].find(j) == dat[i].end())
            dat[i].insert(std::pair<int, int>(j, 1u));
        else
            dat[i][j]++; 
        
        if (!directed)
        {

            if (dat[j].find(i) == dat[j].end())
                dat[j].insert(std::pair<int, int>(i, 1u));
            else
                dat[j][i]++;

        }

        E++;

    }

    N = size;

    return;

}


inline AdjList::AdjList(AdjList && a) :
    dat(std::move(a.dat)),
    directed(a.directed),
    N(a.N),
    E(a.E)
{

}

inline AdjList::AdjList(const AdjList & a) :
    dat(a.dat),
    directed(a.directed),
    N(a.N),
    E(a.E)
{

}

inline AdjList& AdjList::operator=(const AdjList& a)
{
    if (this == &a)
        return *this;

    this->dat = a.dat;
    this->directed = a.directed;
    this->N = a.N;
    this->E = a.E;

    return *this;
}

inline void AdjList::read_edgelist(
    std::string fn,
    int size,
    int skip,
    bool directed
) {

    int i,j;
    std::ifstream filei(fn);

    if (!filei)
        throw std::logic_error("The file " + fn + " was not found.");

    int linenum = 0;
    std::vector< int > source_;
    std::vector< int > target_;

    source_.reserve(1e5);
    target_.reserve(1e5);

    int max_id = size - 1;

    while (!filei.eof())
    {

        if (linenum++ < skip)
            continue;

        filei >> i >> j;

        // Looking for exceptions
        if (filei.bad())
            throw std::logic_error(
                "I/O error while reading the file " +
                fn
            );

        if (filei.fail())
            break;

        if (i > max_id)
            throw std::range_error(
                "The source["+std::to_string(linenum)+"] = " + std::to_string(i) +
                " is above the max_id " + std::to_string(max_id)
                );

        if (j > max_id)
            throw std::range_error(
                "The target["+std::to_string(linenum)+"] = " + std::to_string(j) +
                " is above the max_id " + std::to_string(max_id)
                );

        source_.push_back(i);
        target_.push_back(j);

    }

    // Now using the right constructor
    *this = AdjList(source_, target_, size, directed);

    return;

}

inline std::map<int,int> AdjList::operator()(
    epiworld_fast_uint i
    ) const {

    if (i >= N)
        throw std::range_error(
            "The vertex id " + std::to_string(i) + " is not in the network."
            );

    return dat[i];

}

inline void AdjList::print(epiworld_fast_uint limit) const {


    epiworld_fast_uint counter = 0;
    printf_epiworld("Nodeset:\n");
    int i = -1;
    for (auto & n : dat)
    {

        if (counter++ > limit)
            break;

        printf_epiworld("  % 3i: {", ++i);
        int niter = 0;
        for (auto n_n : n)
            if (++niter < static_cast<int>(n.size()))
            {    
                printf_epiworld("%i, ", static_cast<int>(n_n.first));
            }
            else {
                printf_epiworld("%i}\n", static_cast<int>(n_n.first));
            }
    }

    if (limit < dat.size())
    {
        printf_epiworld(
            "  (... skipping %i records ...)\n",
            static_cast<int>(dat.size() - limit)
            );
    }

}

inline size_t AdjList::vcount() const 
{
    return N;
}

inline size_t AdjList::ecount() const 
{
    return E;
}

inline bool AdjList::is_directed() const {

    if (dat.size() == 0u)
        throw std::logic_error("The edgelist is empty.");
    
    return directed;
    
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/adjlist-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/randgraph.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_RANDGRA
#define EPIWORLD_RANDGRA

template<typename TSeq>
class Model;

template<typename TSeq>
class Agent;

class AdjList;


template<typename TSeq, typename TDat>
inline void rewire_degseq(
    TDat * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    );

template<typename TSeq = int>
inline void rewire_degseq(
    std::vector< Agent<TSeq> > * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    )
{

    #ifdef EPI_DEBUG
    std::vector< int > _degree0(agents->size(), 0);
    for (size_t i = 0u; i < _degree0.size(); ++i)
        _degree0[i] = model->get_agents()[i].get_neighbors().size();
    #endif

    // Identifying individuals with degree > 0
    std::vector< epiworld_fast_uint > non_isolates;
    std::vector< epiworld_double > weights;
    epiworld_double nedges = 0.0;
    
    for (epiworld_fast_uint i = 0u; i < agents->size(); ++i)
    {
        if (agents->operator[](i).get_neighbors().size() > 0u)
        {
            non_isolates.push_back(i);
            epiworld_double wtemp = static_cast<epiworld_double>(
                agents->operator[](i).get_neighbors().size()
                );
            weights.push_back(wtemp);
            nedges += wtemp;
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (epiworld_fast_uint i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    // Only swap if needed
    epiworld_fast_uint N = non_isolates.size();
    epiworld_double prob;
    int nrewires = floor(proportion * nedges);
    while (nrewires-- > 0)
    {

        // Picking egos
        prob = model->runif();
        int id0 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id0 = i;
                break;
            }

        prob = model->runif();
        int id1 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id1 = i;
                break;
            }

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= static_cast<int>(N))
            id1 = 0;

        Agent<TSeq> & p0 = agents->operator[](non_isolates[id0]);
        Agent<TSeq> & p1 = agents->operator[](non_isolates[id1]);

        // Picking alters (relative location in their lists)
        // In this case, these are uniformly distributed within the list
        int id01 = std::floor(p0.get_n_neighbors() * model->runif());
        int id11 = std::floor(p1.get_n_neighbors() * model->runif());

        // When rewiring, we need to flip the individuals from the other
        // end as well, since we are dealing withi an undirected graph
        
        // Finding what neighbour is id0
        model->get_agents()[id0].swap_neighbors(
            model->get_agents()[id1],
            id01,
            id11
            );
        

    }

    #ifdef EPI_DEBUG
    for (size_t _i = 0u; _i < _degree0.size(); ++_i)
    {
        if (_degree0[_i] != static_cast<int>(model->get_agents()[_i].get_n_neighbors()))
            throw std::logic_error("[epi-debug] Degree does not match afted rewire_degseq.");
    }
    #endif

    return;

}

template<typename TSeq>
inline void rewire_degseq(
    AdjList * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    )
{

    // Identifying individuals with degree > 0
    std::vector< epiworld_fast_int > nties(agents->vcount(), 0); 

    #ifdef EPI_DEBUG
    std::vector< int > _degree0(agents->vcount(), 0);
    for (size_t i = 0u; i < _degree0.size(); ++i)
        _degree0[i] = agents->get_dat()[i].size();
    #endif
    
    std::vector< epiworld_fast_uint > non_isolates;
    non_isolates.reserve(nties.size());

    std::vector< epiworld_double > weights;
    weights.reserve(nties.size());

    epiworld_double nedges = 0.0;
    auto & dat = agents->get_dat();

    for (size_t i = 0u; i < dat.size(); ++i)
        nties[i] += dat[i].size();
    
    bool directed = agents->is_directed();
    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (nties[i] > 0)
        {
            non_isolates.push_back(i);
            if (directed)
            {
                weights.push_back( 
                    static_cast<epiworld_double>(nties[i])
                );
                nedges += static_cast<epiworld_double>(nties[i]);
            }
            else {
                weights.push_back( 
                    static_cast<epiworld_double>(nties[i])/2.0
                );
                nedges += static_cast<epiworld_double>(nties[i]) / 2.0;
            }
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (epiworld_fast_uint i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    // Only swap if needed
    epiworld_fast_uint N = non_isolates.size();
    epiworld_double prob;
    int nrewires = floor(proportion * nedges / (
        agents->is_directed() ? 1.0 : 2.0
    ));

    while (nrewires-- > 0)
    {

        // Picking egos
        prob = model->runif();
        int id0 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id0 = i;
                break;
            }

        prob = model->runif();
        int id1 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id1 = i;
                break;
            }

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= static_cast<int>(N))
            id1 = 0;

        std::map<int,int> & p0 = agents->get_dat()[non_isolates[id0]];
        std::map<int,int> & p1 = agents->get_dat()[non_isolates[id1]];

        // Picking alters (relative location in their lists)
        // In this case, these are uniformly distributed within the list
        int id01 = std::floor(p0.size() * model->runif());
        int id11 = std::floor(p1.size() * model->runif());

        // Since it is a map, we need to find the actual ids (positions)
        // are not good enough.
        int count = 0;
        for (auto & n : p0)
            if (count++ == id01)
                id01 = n.first;

        count = 0;
        for (auto & n : p1)
            if (count++ == id11)
                id11 = n.first;

        // When rewiring, we need to flip the individuals from the other
        // end as well, since we are dealing withi an undirected graph
        
        // Finding what neighbour is id0
        if (!agents->is_directed())
        {

            std::map<int,int> & p01 = agents->get_dat()[id01];
            std::map<int,int> & p11 = agents->get_dat()[id11];

            std::swap(p01[id0], p11[id1]);
            
        }

        // Moving alter first
        std::swap(p0[id01], p1[id11]);

    }

    #ifdef EPI_DEBUG
    for (size_t _i = 0u; _i < _degree0.size(); ++_i)
    {
        if (_degree0[_i] != static_cast<int>(agents->get_dat()[_i].size()))
            throw std::logic_error(
                "[epi-debug] Degree does not match afted rewire_degseq. " +
                std::string("Expected: ") + 
                std::to_string(_degree0[_i]) + 
                std::string(", observed: ") +
                std::to_string(agents->get_dat()[_i].size())
                );
    }
    #endif


    return;

}

template<typename TSeq>
inline AdjList rgraph_bernoulli(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    std::vector< int > source;
    std::vector< int > target;

    // Checking the density (how many)
    std::binomial_distribution<> d(
        n * (n - 1.0) / (directed ? 1.0 : 2.0),
        p
    );

    epiworld_fast_uint m = d(model.get_rand_endgine());

    source.resize(m);
    target.resize(m);

    epiworld_fast_uint a,b;
    for (epiworld_fast_uint i = 0u; i < m; ++i)
    {
        a = floor(model.runif() * n);

        if (!directed)
            b = floor(model.runif() * a);
        else
        {
            b = floor(model.runif() * n);
            if (b == a)
                b++;
            
            if (b >= n)
                b = 0u;
        }

        source[i] = static_cast<int>(a);
        target[i] = static_cast<int>(b);

    }

    AdjList al(source, target, static_cast<int>(n), directed);

    return al;
    
}

template<typename TSeq>
inline AdjList rgraph_bernoulli2(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    std::vector< int > source;
    std::vector< int > target;

    // Checking the density (how many)
    std::binomial_distribution<> d(
        n * (n - 1.0) / (directed ? 1.0 : 2.0),
        p
    );

    // Need to compensate for the possible number of diagonal
    // elements sampled. If n * n, then each diag element has
    // 1/(n^2) chance of sampling

    epiworld_fast_uint m = d(model.get_rand_endgine());

    source.resize(m);
    target.resize(m);

    double n2 = static_cast<double>(n * n);

    int loc,row,col;
    for (epiworld_fast_uint i = 0u; i < m; ++i)
    {
        loc = floor(model.runif() * n2);
        col = floor(static_cast<double>(loc)/static_cast<double>(n));
        row = loc - row * n;

        // Undirected needs to swap
        if (!directed && (col > row))
            std::swap(col, row);

        source[i] = row;
        target[i] = col;

    }

    AdjList al(source, target, static_cast<int>(n), directed);

    return al;
    
}

inline AdjList rgraph_ring_lattice(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    bool directed = false
) {

    if ((n - 1u) < k)
        throw std::logic_error("k can be at most n - 1.");

    std::vector< int > source;
    std::vector< int > target;

    if (!directed)
        if (k > 1u) k = static_cast< size_t >(floor(k / 2.0));

    for (size_t i = 0; i < n; ++i)
    {

        for (size_t j = 1u; j <= k; ++j)
        {

            // Next neighbor
            size_t l = i + j;
            if (l >= n) l = l - n;

            source.push_back(i);
            target.push_back(l);

        }

    }

    return AdjList(source, target, n, directed);

}

/**
 * @brief Smallworld network (Watts-Strogatz)
 * 
 * @tparam TSeq 
 * @param n 
 * @param k 
 * @param p 
 * @param directed 
 * @param model 
 * @return AdjList 
 */
template<typename TSeq>
inline AdjList rgraph_smallworld(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    // Creating the ring lattice
    AdjList ring = rgraph_ring_lattice(n,k,directed);
    
    // Rewiring and returning
    if (k > 0u)
        rewire_degseq(&ring, &model, p);
        
    return ring;

}

/**
 * @brief Generates a blocked network
 * 
 * Since block sizes and number of connections between blocks are fixed,
 * this routine is fully deterministic.
 * 
 * @tparam TSeq 
 * @param n Size of the network
 * @param blocksize Size of the block.
 * @param ncons Number of connections between blocks
 * @param model A model
 * @return AdjList 
 */
template<typename TSeq>
inline AdjList rgraph_blocked(
    epiworld_fast_uint n,
    epiworld_fast_uint blocksize,
    epiworld_fast_uint ncons,
    Model<TSeq>&
) {

    std::vector< int > source_;
    std::vector< int > target_;

    size_t i = 0u;
    size_t cum_node_count = 0u;
    while (i < n)
    {

        for (size_t j = 0; j < blocksize; ++j)
        {

            for (size_t k = 0; k < j; ++k)
            {
                // No loops
                if (k == j)
                    continue;

                // Exists the loop in case there are no more 
                // nodes available
                if ((i + k) >= n)
                    break;

                source_.push_back(static_cast<int>(j + i));
                target_.push_back(static_cast<int>(k + i));
            }

            // No more nodes left to build connections
            if (++cum_node_count >= n)
                break;
            
        }

        // Connections between this and the previou sone
        if (i != 0)
        {

            size_t max_cons = std::min(ncons, n - cum_node_count);

            // Generating the connections
            for (size_t j = 0u; j < max_cons; ++j)
            {

                source_.push_back(static_cast<int>(i + j - blocksize));
                target_.push_back(static_cast<int>(i + j));

            }
        }

        i += blocksize;
        
    }
        
    return AdjList(source_, target_, n, false);

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/randgraph.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/queue-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_QUEUE_BONES_HPP
#define EPIWORLD_QUEUE_BONES_HPP

/**
 * @brief Controls which agents are verified at each step
 * 
 * @details The idea is that only agents who are either in
 * an infected state or have an infected neighbor should be
 * checked. Otherwise it makes no sense (no chance to recover
 * or capture the disease).
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Queue
{
    friend class Model<TSeq>;

private:

    /**
     * @brief Count of ego's neighbors in queue (including ego)
     */
    std::vector< epiworld_fast_int > active;
    Model<TSeq> * model = nullptr;
    int n_in_queue = 0;

    // Auxiliary variable that checks how many steps
    // left are there
    // int n_steps_left;
    // bool queuing_started   = false;

public:

    void operator+=(Agent<TSeq> * p);
    void operator-=(Agent<TSeq> * p);
    epiworld_fast_int & operator[](epiworld_fast_uint i);

    // void initialize(Model<TSeq> * m, Agent<TSeq> * p);
    void reset();

    bool operator==(const Queue<TSeq> & other) const;
    bool operator!=(const Queue<TSeq> & other) const {return !operator==(other);};

};

template<typename TSeq>
inline void Queue<TSeq>::operator+=(Agent<TSeq> * p)
{

    if (++active[p->id] == 1)
        n_in_queue++;

    for (auto n : p->neighbors)
    {

        if (++active[n] == 1)
            n_in_queue++;

    }

}

template<typename TSeq>
inline void Queue<TSeq>::operator-=(Agent<TSeq> * p)
{

    if (--active[p->id] == 0)
        n_in_queue--;

    for (auto n : p->neighbors)
    {
        if (--active[n] == 0)
            n_in_queue--;
    }

}

template<typename TSeq>
inline epiworld_fast_int & Queue<TSeq>::operator[](epiworld_fast_uint i)
{
    return active[i];
}

template<typename TSeq>
inline void Queue<TSeq>::reset()
{

    if (n_in_queue)
    {

        for (auto & q : this->active)
            q = 0;

        n_in_queue = 0;
        
    }

    active.resize(model->size(), 0);

}

template<typename TSeq>
inline bool Queue<TSeq>::operator==(const Queue<TSeq> & other) const 
{
    if (active.size() != other.active.size())
        return false;

    for (size_t i = 0u; i < active.size(); ++i)
    {
        if (active[i] != other.active[i])
            return false;
    }

    return true;
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/queue-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/model-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODEL_HPP
#define EPIWORLD_MODEL_HPP

template<typename TSeq>
class Agent;

template<typename TSeq>
class AgentsSample;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Tool;

class AdjList;

template<typename TSeq>
class DataBase;

template<typename TSeq>
class Queue;

template<typename TSeq>
struct Action;

template<typename TSeq>
inline epiworld_double susceptibility_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double transmission_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double recovery_enhancer_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double death_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );

template<typename TSeq>
inline std::function<void(size_t,Model<TSeq>*)> make_save_run(
    std::string fmt = "%03lu-episimulation.csv",
    bool total_hist = true,
    bool variant_info = false,
    bool variant_hist = false,
    bool tool_info = false,
    bool tool_hist = false,
    bool transmission = false,
    bool transition = false,
    bool reproductive = false
    );

// template<typename TSeq>
// class VirusPtr;

// template<typename TSeq>
// class ToolPtr;

/**
 * @brief Core class of epiworld.
 * 
 * The model class provides the wrapper that puts together `Agent`, `Virus`, and
 * `Tools`.
 * 
 * @tparam TSeq Type of sequence. In principle, users can build models in which
 * virus and human sequence is represented as numeric vectors (if needed.)
 */
template<typename TSeq>
class Model {
    friend class Agent<TSeq>;
    friend class AgentsSample<TSeq>;
    friend class DataBase<TSeq>;
    friend class Queue<TSeq>;
protected:

    std::string name = ""; ///< Name of the model

    DataBase<TSeq> db = DataBase<TSeq>(*this);

    std::vector< Agent<TSeq> > population = {};

    bool using_backup = true;
    std::vector< Agent<TSeq> > population_backup = {};


    /**
     * @name Auxiliary variables for AgentsSample<TSeq> iterators
     * 
     * @details These variables+objects are used by the AgentsSample<TSeq>
     * class for building efficient iterators over agents. The idea is to
     * reduce the memory allocation, so only during the first call of
     * AgentsSample<TSeq>::AgentsSample(Model<TSeq>) these vectors are allocated.
     */
    ///@{
    std::vector< Agent<TSeq> * > sampled_population;
    size_t sampled_population_n = 0u;
    std::vector< size_t > population_left;
    size_t population_left_n = 0u;
    ///@}

    /**
     * @name Agents features
     * 
     * @details Optionally, a model can include an external data source
     * pointing to agents information. The data can then be access through
     * the `Agent::operator()` method.
     * 
     */
    ///@{
    double * agents_data = nullptr;
    size_t agents_data_ncols = 0u;
    ///@}

    bool directed = false;
    
    std::vector< VirusPtr<TSeq> > viruses = {};
    std::vector< epiworld_double > prevalence_virus = {}; ///< Initial prevalence_virus of each virus
    std::vector< bool > prevalence_virus_as_proportion = {};
    std::vector< VirusToAgentFun<TSeq> > viruses_dist_funs = {};
    
    std::vector< ToolPtr<TSeq> > tools = {};
    std::vector< epiworld_double > prevalence_tool = {};
    std::vector< bool > prevalence_tool_as_proportion = {};
    std::vector< ToolToAgentFun<TSeq> > tools_dist_funs = {};

    std::vector< Entity<TSeq> > entities = {}; 
    std::vector< Entity<TSeq> > entities_backup = {};

    std::mt19937 engine;
    
    std::uniform_real_distribution<> runifd      = std::uniform_real_distribution<> (0.0, 1.0);
    std::normal_distribution<>       rnormd      = std::normal_distribution<>(0.0);
    std::gamma_distribution<>        rgammad     = std::gamma_distribution<>();
    std::lognormal_distribution<>    rlognormald = std::lognormal_distribution<>();
    std::exponential_distribution<>  rexpd       = std::exponential_distribution<>();

    std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> rewire_fun;
    epiworld_double rewire_prop = 0.0;
        
    std::map<std::string, epiworld_double > parameters;
    epiworld_fast_uint ndays = 0;
    Progress pb;

    std::vector< UpdateFun<TSeq> >    status_fun = {};
    std::vector< std::string >        states_labels = {};
    epiworld_fast_uint nstatus = 0u;
    
    bool verbose     = true;
    int current_date = 0;

    void dist_tools();
    void dist_virus();
    // void dist_entities();

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed = 
        std::chrono::duration<epiworld_double,std::micro>::zero();
    epiworld_fast_uint n_replicates = 0u;
    void chrono_start();
    void chrono_end();

    std::vector<std::function<void(Model<TSeq>*)>> global_action_functions;
    std::vector< int > global_action_dates;

    Queue<TSeq> queue;
    bool use_queuing   = true;

    /**
     * @brief Variables used to keep track of the actions
     * to be made regarding viruses.
     */
    std::vector< Action<TSeq> > actions = {};
    epiworld_fast_uint nactions = 0u;

    /**
     * @brief Construct a new Action object
     * 
     * @param agent_ Agent over which the action will be called
     * @param virus_ Virus pointer included in the action
     * @param tool_ Tool pointer included in the action
     * @param entity_ Entity pointer included in the action
     * @param new_state_ New state of the agent
     * @param call_ Function the action will call
     * @param queue_ Change in the queue
     * @param idx_agent_ Location of agent in object.
     * @param idx_object_ Location of object in agent.
     */
    void actions_add(
        Agent<TSeq> * agent_,
        VirusPtr<TSeq> virus_,
        ToolPtr<TSeq> tool_,
        Entity<TSeq> * entity_,
        epiworld_fast_uint new_state_,
        epiworld_fast_int queue_,
        ActionFun<TSeq> call_,
        int idx_agent_,
        int idx_object_
        );

    /**
     * @brief Executes the stored action
     * 
     * @param model_ Model over which it will be executed.
     */
    void actions_run();

    /**
     * @name Tool Mixers
     * 
     * These functions combine the effects tools have to deliver
     * a single effect. For example, wearing a mask, been vaccinated,
     * and the immune system combine together to jointly reduce
     * the susceptibility for a given virus.
     * 
     */
    MixerFun<TSeq> susceptibility_reduction_mixer = susceptibility_reduction_mixer_default<TSeq>;
    MixerFun<TSeq> transmission_reduction_mixer = transmission_reduction_mixer_default<TSeq>;
    MixerFun<TSeq> recovery_enhancer_mixer = recovery_enhancer_mixer_default<TSeq>;
    MixerFun<TSeq> death_reduction_mixer = death_reduction_mixer_default<TSeq>;

    /**
     * @brief Advanced usage: Makes a copy of data and returns it as undeleted pointer
     * 
     * @param copy 
     */
    virtual Model<TSeq> * clone_ptr();

public:

    std::vector<epiworld_double> array_double_tmp;
    std::vector<Virus<TSeq> * > array_virus_tmp;

    Model();
    Model(const Model<TSeq> & m);
    Model(Model<TSeq> & m) = delete;
    Model(Model<TSeq> && m);
    Model<TSeq> & operator=(const Model<TSeq> & m);

    virtual ~Model() {};

    void clone_population(
        std::vector< Agent<TSeq> > & other_population,
        std::vector< Entity<TSeq> > & other_entities,
        Model<TSeq> * other_model,
        bool & other_directed
    ) const ;

    void clone_population(
        const Model<TSeq> & other_model
    );

    /**
     * @name Set the backup object
     * @details `backup` can be used to restore the entire object
     * after a run. This can be useful if the user wishes to have
     * individuals start with the same network from the beginning.
     * 
     */
    ///@{
    void set_backup();
    // void restore_backup();
    ///@}

    DataBase<TSeq> & get_db();
    epiworld_double & operator()(std::string pname);

    size_t size() const;

    /**
     * @name Random number generation
     * 
     * @param eng Random number generator
     * @param s Seed
     */
    ///@{
    void set_rand_engine(std::mt19937 & eng);
    std::mt19937 & get_rand_endgine();
    void seed(size_t s);
    void set_rand_norm(epiworld_double mean, epiworld_double sd);
    void set_rand_unif(epiworld_double a, epiworld_double b);
    void set_rand_exp(epiworld_double lambda);
    void set_rand_gamma(epiworld_double alpha, epiworld_double beta);
    void set_rand_lognormal(epiworld_double mean, epiworld_double shape);
    epiworld_double runif();
    epiworld_double runif(epiworld_double a, epiworld_double b);
    epiworld_double rnorm();
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma();
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double rexp();
    epiworld_double rexp(epiworld_double lambda);
    epiworld_double rlognormal();
    epiworld_double rlognormal(epiworld_double mean, epiworld_double shape);
    ///@}

    /**
     * @name Add Virus/Tool to the model
     * 
     * This is done before the model has been initialized.
     * 
     * @param v Virus to be added
     * @param t Tool to be added
     * @param preval Initial prevalence (initial state.) It can be
     * specified as a proportion (between zero and one,) or an integer
     * indicating number of individuals.
     */
    ///@{
    void add_virus(Virus<TSeq> v, epiworld_double preval);
    void add_virus_n(Virus<TSeq> v, epiworld_fast_uint preval);
    void add_virus_fun(Virus<TSeq> v, VirusToAgentFun<TSeq> fun);
    void add_tool(Tool<TSeq> t, epiworld_double preval);
    void add_tool_n(Tool<TSeq> t, epiworld_fast_uint preval);
    void add_tool_fun(Tool<TSeq> t, ToolToAgentFun<TSeq> fun);
    void add_entity(Entity<TSeq> e);
    void rm_virus(size_t virus_pos);
    void rm_tool(size_t tool_pos);
    void rm_entity(size_t entity_pos);
    ///@}

    /**
     * @brief Associate agents-entities from a file
     * 
     * The structure of the file should be two columns separated by 
     * space. The first column indexing between 0 and nagents-1, and the
     * second column between 0 and nentities - 1.
     * 
     * @param fn Path to the file.
     * @param skip How many rows to skip.
     */
    void load_agents_entities_ties(std::string fn, int skip);

    /**
     * @name Accessing population of the model
     * 
     * @param fn std::string Filename of the edgelist file.
     * @param skip int Number of lines to skip in `fn`.
     * @param directed bool Whether the graph is directed or not.
     * @param size Size of the network.
     * @param al AdjList to read into the model.
     */
    ///@{
    void agents_from_adjlist(
        std::string fn,
        int size,
        int skip = 0,
        bool directed = false
        );

    void agents_from_edgelist(
        const std::vector< int > & source,
        const std::vector< int > & target,
        int size,
        bool directed
    );

    void agents_from_adjlist(AdjList al);

    bool is_directed() const;

    std::vector< Agent<TSeq> > & get_agents();

    std::vector< Entity<TSeq> > & get_entities();

    void agents_smallworld(
        epiworld_fast_uint n = 1000,
        epiworld_fast_uint k = 5,
        bool d = false,
        epiworld_double p = .01
        );
    void agents_empty_graph(epiworld_fast_uint n = 1000);
    ///@}

    /**
     * @name Functions to run the model
     * 
     * @param seed Seed to be used for Pseudo-RNG.
     * @param ndays Number of days (steps) of the simulation.
     * @param fun In the case of `run_multiple`, a function that is called
     * after each experiment.
     * 
     */
    ///@{
    void update_state();
    void mutate_variant();
    void next();
    virtual void run(
        epiworld_fast_uint ndays,
        int seed = -1
    ); ///< Runs the simulation (after initialization)
    void run_multiple( ///< Multiple runs of the simulation
        epiworld_fast_uint ndays,
        epiworld_fast_uint nexperiments,
        int seed_ = -1,
        std::function<void(size_t,Model<TSeq>*)> fun = make_save_run<TSeq>(),
        bool reset = true,
        bool verbose = true,
        int nthreads = 1
        );
    ///@}

    size_t get_n_variants() const;
    size_t get_n_tools() const;
    epiworld_fast_uint get_ndays() const;
    epiworld_fast_uint get_n_replicates() const;
    void set_ndays(epiworld_fast_uint ndays);
    bool get_verbose() const;
    void verbose_off();
    void verbose_on();
    int today() const; ///< The current time of the model

    /**
     * @name Rewire the network preserving the degree sequence.
     *
     * @details This implementation assumes an undirected network,
     * thus if {(i,j), (k,l)} -> {(i,l), (k,j)}, the reciprocal
     * is also true, i.e., {(j,i), (l,k)} -> {(j,k), (l,i)}.
     * 
     * @param proportion Proportion of ties to be rewired.
     * 
     * @result A rewired version of the network.
     */
    ///@{
    void set_rewire_fun(std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> fun);
    void set_rewire_prop(epiworld_double prop);
    epiworld_double get_rewire_prop() const;
    void rewire();
    ///@}

    /**
     * @brief Wrapper of `DataBase::write_data`
     * 
     * @param fn_variant_info Filename. Information about the variant.
     * @param fn_variant_hist Filename. History of the variant.
     * @param fn_tool_info Filename. Information about the tool.
     * @param fn_tool_hist Filename. History of the tool.
     * @param fn_total_hist   Filename. Aggregated history (state)
     * @param fn_transmission Filename. Transmission history.
     * @param fn_transition   Filename. Markov transition history.
     * @param fn_reproductive_number Filename. Case by case reproductive number
     */
    void write_data(
        std::string fn_variant_info,
        std::string fn_variant_hist,
        std::string fn_tool_info,
        std::string fn_tool_hist,
        std::string fn_total_hist,
        std::string fn_transmission,
        std::string fn_transition,
        std::string fn_reproductive_number,
        std::string fn_generation_time
        ) const;

    /**
     * @name Export the network data in edgelist form
     * 
     * @param fn std::string. File name.
     * @param source Integer vector
     * @param target Integer vector
     * 
     * @details When passing the source and target, the function will
     * write the edgelist on those.
     */
    ///@{
    void write_edgelist(
        std::string fn
        ) const;

    void write_edgelist(
        std::vector< epiworld_fast_uint > & source,
        std::vector< epiworld_fast_uint > & target
        ) const;
    ///@}

    std::map<std::string, epiworld_double> & params();

    /**
     * @brief Reset the model
     * 
     * @details Resetting the model will:
     * - clear the database
     * - restore the population (if `set_backup()` was called before)
     * - re-distribute tools
     * - re-distribute viruses
     * - set the date to 0
     * 
     */
    void reset();
    void print(bool lite = false) const;

    Model<TSeq> && clone() const;

    /**
     * @name Manage state (states) in the model
     * 
     * @details
     * 
     * The functions `get_state` return the current values for the 
     * states included in the model.
     * 
     * @param lab `std::string` Name of the state.
     * 
     * @return `add_state*` returns nothing.
     * @return `get_state_*` returns a vector of pairs with the 
     * states and their labels.
     */
    ///@{
    void add_state(std::string lab, UpdateFun<TSeq> fun = nullptr);
    const std::vector< std::string > & get_state() const;
    const std::vector< UpdateFun<TSeq> > & get_state_fun() const;
    void print_state_codes() const;
    ///@}

    /**
     * @name Setting and accessing parameters from the model
     * 
     * @details Tools can incorporate parameters included in the model.
     * Internally, parameters in the tool are stored as pointers to
     * an std::map<> of parameters in the model. Using the `epiworld_fast_uint`
     * method directly fetches the parameters in the order these were
     * added to the tool. Accessing parameters via the `std::string` method
     * involves searching the parameter directly in the std::map<> member
     * of the model (so it is not recommended.)
     * 
     * The `par()` function members are aliases for `get_param()`.
     * 
     * In the case of the function `read_params`, users can pass a file
     * listing parameters to be included in the model. Each line in the
     * file should have the following structure:
     * 
     * ```
     * [name of parameter 1]: [value in double]
     * [name of parameter 2]: [value in double]
     * ...
     * ```
     * 
     * The only condition for parameter names is that these do not include
     * a colon.
     * 
     * 
     * @param initial_val 
     * @param pname Name of the parameter to add or to fetch
     * @param fn Path to the file containing parameters
     * @return The current value of the parameter
     * in the model.
     * 
     */
    ///@{
    epiworld_double add_param(epiworld_double initial_val, std::string pname);
    void read_params(std::string fn);
    epiworld_double get_param(epiworld_fast_uint k);
    epiworld_double get_param(std::string pname);
    epiworld_double par(epiworld_fast_uint k);
    epiworld_double par(std::string pname);
    ///@}

    void get_elapsed(
        std::string unit = "auto",
        epiworld_double * last_elapsed = nullptr,
        epiworld_double * total_elapsed = nullptr,
        std::string * unit_abbr = nullptr,
        bool print = true
    ) const;

    /**
     * @name Set the user data object
     * 
     * @param names string vector with the names of the variables.
     */
    ///[@
    void set_user_data(std::vector< std::string > names);
    void add_user_data(epiworld_fast_uint j, epiworld_double x);
    void add_user_data(std::vector< epiworld_double > x);
    UserData<TSeq> & get_user_data();
    ///@}

    /**
     * @brief Set a global action
     * 
     * @param fun A function to be called on the prescribed dates
     * @param date Integer indicating when the function is called (see details)
     * 
     * @details When date is less than zero, then the function is called
     * at the end of every day. Otherwise, the function will be called only
     * at the end of the indicated date.
     */
    void add_global_action(
        std::function<void(Model<TSeq>*)> fun,
        int date = -99
        );

    void run_global_actions();

    void clear_state_set();

    /**
     * @name Queuing system
     * @details When queueing is on, the model will keep track of which agents
     * are either in risk of exposure or exposed. This then is used at each 
     * step to act only on the aforementioned agents.
     * 
     */
    ////@{
    void queuing_on(); ///< Activates the queuing system (default.)
    void queuing_off(); ///< Deactivates the queuing system.
    bool is_queuing_on() const; ///< Query if the queuing system is on.
    Queue<TSeq> & get_queue(); ///< Retrieve the `Queue` object.
    ///@}

    /**
     * @name Get the susceptibility reduction object
     * 
     * @param v 
     * @return epiworld_double 
     */
    ///@{
    void set_susceptibility_reduction_mixer(MixerFun<TSeq> fun);
    void set_transmission_reduction_mixer(MixerFun<TSeq> fun);
    void set_recovery_enhancer_mixer(MixerFun<TSeq> fun);
    void set_death_reduction_mixer(MixerFun<TSeq> fun);
    ///@}

    const std::vector< VirusPtr<TSeq> > & get_viruses() const;
    const std::vector< ToolPtr<TSeq> > & get_tools() const;

    /**
     * @brief Set the agents data object
     * 
     * @details The data should be an array with the data stored in a
     * column major order, i.e., by column.
     * 
     * @param data_ Pointer to the first element of an array of size
     * `size() * ncols_`.
     * @param ncols_ Number of features included in the data.
     * 
     */
    void set_agents_data(double * data_, size_t ncols_);
    double * get_agents_data();
    size_t get_agents_data_ncols();

    /**
     * @brief Set the name object
     * 
     * @param name 
     */
    void set_name(std::string name);
    std::string get_name() const;

    bool operator==(const Model<TSeq> & other) const;
    bool operator!=(const Model<TSeq> & other) const {return !operator==(other);};

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/model-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/model-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODEL_MEAT_HPP
#define EPIWORLD_MODEL_MEAT_HPP



/**
 * @brief Function factory for saving model runs
 * 
 * @details This function is the default behavior of the `run_multiple`
 * member of `Model<TSeq>`. By default only the total history (
 * case counts by unit of time.)
 * 
 * @tparam TSeq 
 * @param fmt 
 * @param total_hist 
 * @param variant_info 
 * @param variant_hist 
 * @param tool_info 
 * @param tool_hist 
 * @param transmission 
 * @param transition 
 * @return std::function<void(size_t,Model<TSeq>*)> 
 */
template<typename TSeq = int>
inline std::function<void(size_t,Model<TSeq>*)> make_save_run(
    std::string fmt,
    bool total_hist,
    bool variant_info,
    bool variant_hist,
    bool tool_info,
    bool tool_hist,
    bool transmission,
    bool transition,
    bool reproductive,
    bool generation
    )
{

    // Counting number of %
    int n_fmt = 0;
    for (auto & f : fmt)
        if (f == '%')
            n_fmt++;

    if (n_fmt != 1)
        throw std::logic_error("The -fmt- argument must have only one \"%\" symbol.");

    // Listting things to save
    std::vector< bool > what_to_save = {
        variant_info,
        variant_hist,
        tool_info,
        tool_hist,
        total_hist,
        transmission,
        transition,
        reproductive,
        generation
    };

    std::function<void(size_t,Model<TSeq>*)> saver = [fmt,what_to_save](
        size_t niter, Model<TSeq> * m
    ) -> void {

        std::string variant_info = "";
        std::string variant_hist = "";
        std::string tool_info = "";
        std::string tool_hist = "";
        std::string total_hist = "";
        std::string transmission = "";
        std::string transition = "";
        std::string reproductive = "";
        std::string generation = "";

        char buff[128];
        if (what_to_save[0u])
        {
            variant_info = fmt + std::string("_variant_info.csv");
            snprintf(buff, sizeof(buff), variant_info.c_str(), niter);
            variant_info = buff;
        } 
        if (what_to_save[1u])
        {
            variant_hist = fmt + std::string("_variant_hist.csv");
            snprintf(buff, sizeof(buff), variant_hist.c_str(), niter);
            variant_hist = buff;
        } 
        if (what_to_save[2u])
        {
            tool_info = fmt + std::string("_tool_info.csv");
            snprintf(buff, sizeof(buff), tool_info.c_str(), niter);
            tool_info = buff;
        } 
        if (what_to_save[3u])
        {
            tool_hist = fmt + std::string("_tool_hist.csv");
            snprintf(buff, sizeof(buff), tool_hist.c_str(), niter);
            tool_hist = buff;
        } 
        if (what_to_save[4u])
        {
            total_hist = fmt + std::string("_total_hist.csv");
            snprintf(buff, sizeof(buff), total_hist.c_str(), niter);
            total_hist = buff;
        } 
        if (what_to_save[5u])
        {
            transmission = fmt + std::string("_transmission.csv");
            snprintf(buff, sizeof(buff), transmission.c_str(), niter);
            transmission = buff;
        } 
        if (what_to_save[6u])
        {
            transition = fmt + std::string("_transition.csv");
            snprintf(buff, sizeof(buff), transition.c_str(), niter);
            transition = buff;
        } 
        if (what_to_save[7u])
        {

            reproductive = fmt + std::string("_reproductive.csv");
            snprintf(buff, sizeof(buff), reproductive.c_str(), niter);
            reproductive = buff;

        }
        if (what_to_save[8u])
        {

            generation = fmt + std::string("_generation.csv");
            snprintf(buff, sizeof(buff), generation.c_str(), niter);
            generation = buff;

        }
        
    
        m->write_data(
            variant_info,
            variant_hist,
            tool_info,
            tool_hist,
            total_hist,
            transmission,
            transition,
            reproductive,
            generation
        );

    };

    return saver;
}

template<typename TSeq>
inline void Model<TSeq>::actions_add(
    Agent<TSeq> * agent_,
    VirusPtr<TSeq> virus_,
    ToolPtr<TSeq> tool_,
    Entity<TSeq> * entity_,
    epiworld_fast_uint new_state_,
    epiworld_fast_int queue_,
    ActionFun<TSeq> call_,
    int idx_agent_,
    int idx_object_
) {
    
    ++nactions;

    #ifdef EPI_DEBUG
    if (nactions == 0)
        throw std::logic_error("Actions cannot be zero!!");

    if ((virus_ != nullptr) && idx_agent_ >= 0)
    {
        if (idx_agent_ >= static_cast<int>(virus_->get_agent()->get_n_viruses()))
            throw std::logic_error(
                "The virus to add is out of range in the host agent."
                );
    }
    #endif

    if (nactions > actions.size())
    {

        actions.push_back(
            Action<TSeq>(
                agent_, virus_, tool_, entity_, new_state_, queue_, call_,
                idx_agent_, idx_object_
            ));

    }
    else 
    {

        Action<TSeq> & A = actions.at(nactions - 1u);

        A.agent      = agent_;
        A.virus      = virus_;
        A.tool       = tool_;
        A.entity     = entity_;
        A.new_state = new_state_;
        A.queue      = queue_;
        A.call       = call_;
        A.idx_agent  = idx_agent_;
        A.idx_object = idx_object_;

    }

    return;

}

template<typename TSeq>
inline void Model<TSeq>::actions_run()
{
    // Making the call
    while (nactions != 0u)
    {

        Action<TSeq>   a = actions[--nactions];
        Agent<TSeq> * p  = a.agent;

        // Applying function
        if (a.call)
        {
            a.call(a, this);
        }

        // Updating state
        if (static_cast<epiworld_fast_int>(p->state) != a.new_state)
        {

            if (a.new_state >= static_cast<epiworld_fast_int>(nstatus))
                throw std::range_error(
                    "The proposed state " + std::to_string(a.new_state) + " is out of range. " +
                    "The model currently has " + std::to_string(nstatus - 1) + " states.");

            // Figuring out if we need to undo a change
            // If the agent has made a change in the state recently, then we
            // need to undo the accounting, e.g., if A->B was made, we need to
            // undo it and set B->A so that the daily accounting is right.
            if (p->status_last_changed == today())
            {

                // Updating accounting
                db.update_state(p->status_prev, p->state, true); // Undoing
                db.update_state(p->status_prev, a.new_state);

                for (size_t v = 0u; v < p->n_viruses; ++v)
                {
                    db.update_virus(p->viruses[v]->id, p->state, p->status_prev); // Undoing
                    db.update_virus(p->viruses[v]->id, p->status_prev, a.new_state);
                }

                for (size_t t = 0u; t < p->n_tools; ++t)
                {
                    db.update_tool(p->tools[t]->id, p->state, p->status_prev); // Undoing
                    db.update_tool(p->tools[t]->id, p->status_prev, a.new_state);
                }

                // Changing to the new state, we won't update the
                // previous state in case we need to undo the change
                p->state = a.new_state;

            } else {

                // Updating accounting
                db.update_state(p->state, a.new_state);

                for (size_t v = 0u; v < p->n_viruses; ++v)
                    db.update_virus(p->viruses[v]->id, p->state, a.new_state);

                for (size_t t = 0u; t < p->n_tools; ++t)
                    db.update_tool(p->tools[t]->id, p->state, a.new_state);


                // Saving the last state and setting the new one
                p->status_prev = p->state;
                p->state      = a.new_state;

                // It used to be a day before, but we still
                p->status_last_changed = today();

            }
            
        }

        #ifdef EPI_DEBUG
        if (static_cast<int>(p->state) >= static_cast<int>(nstatus))
                throw std::range_error(
                    "The new state " + std::to_string(p->state) + " is out of range. " +
                    "The model currently has " + std::to_string(nstatus - 1) + " states.");
        #endif

        // Updating queue
        if (use_queuing)
        {

            if (a.queue == QueueValues::Everyone)
                queue += p;
            else if (a.queue == -QueueValues::Everyone)
                queue -= p;
            else if (a.queue == QueueValues::OnlySelf)
                queue[p->get_id()]++;
            else if (a.queue == -QueueValues::OnlySelf)
                queue[p->get_id()]--;
            else if (a.queue != QueueValues::NoOne)
                throw std::logic_error(
                    "The proposed queue change is not valid. Queue values can be {-2, -1, 0, 1, 2}."
                    );
                    
        }

    }

    return;
    
}

/**
 * @name Default function for combining susceptibility_reduction levels
 * 
 * @tparam TSeq 
 * @param pt 
 * @return epiworld_double 
 */
///@{
template<typename TSeq>
inline epiworld_double susceptibility_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq> * m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_susceptibility_reduction(v, m));

    return 1.0 - total;
    
}

template<typename TSeq>
inline epiworld_double transmission_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_transmission_reduction(v, m));

    return (1.0 - total);
    
}

template<typename TSeq>
inline epiworld_double recovery_enhancer_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_recovery_enhancer(v, m));

    return 1.0 - total;
    
}

template<typename TSeq>
inline epiworld_double death_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
) {

    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
    {
        total *= (1.0 - tool->get_death_reduction(v, m));
    } 

    return 1.0 - total;
    
}
///@}

template<typename TSeq>
inline Model<TSeq> * Model<TSeq>::clone_ptr()
{
    Model<TSeq> * ptr = new Model<TSeq>(*dynamic_cast<const Model<TSeq>*>(this));

    #ifdef EPI_DEBUG
    if (*this != *ptr)
        throw std::logic_error("Model::clone_ptr The copies of the model don't match.");
    #endif

    return ptr;
}

template<typename TSeq>
inline Model<TSeq>::Model()
{
    db.model = this;
    db.user_data = this;
    if (use_queuing)
        queue.model = this;
}

template<typename TSeq>
inline Model<TSeq>::Model(const Model<TSeq> & model) :
    name(model.name),
    db(model.db),
    population(model.population),
    population_backup(model.population_backup),
    directed(model.directed),
    viruses(model.viruses),
    prevalence_virus(model.prevalence_virus),
    prevalence_virus_as_proportion(model.prevalence_virus_as_proportion),
    viruses_dist_funs(model.viruses_dist_funs),
    tools(model.tools),
    prevalence_tool(model.prevalence_tool),
    prevalence_tool_as_proportion(model.prevalence_tool_as_proportion),
    tools_dist_funs(model.tools_dist_funs),
    entities(model.entities),
    entities_backup(model.entities_backup),
    // prevalence_entity(model.prevalence_entity),
    // prevalence_entity_as_proportion(model.prevalence_entity_as_proportion),
    // entities_dist_funs(model.entities_dist_funs),
    rewire_fun(model.rewire_fun),
    rewire_prop(model.rewire_prop),
    parameters(model.parameters),
    ndays(model.ndays),
    pb(model.pb),
    status_fun(model.status_fun),
    states_labels(model.states_labels),
    nstatus(model.nstatus),
    verbose(model.verbose),
    current_date(model.current_date),
    global_action_functions(model.global_action_functions),
    global_action_dates(model.global_action_dates),
    queue(model.queue),
    use_queuing(model.use_queuing),
    array_double_tmp(model.array_double_tmp.size()),
    array_virus_tmp(model.array_virus_tmp.size())
{


    // Removing old neighbors
    for (auto & p : population)
        p.model = this;

    if (population_backup.size() != 0u)
        for (auto & p : population_backup)
            p.model = this;

    for (auto & e : entities)
        e.model = this;

    if (entities_backup.size() != 0u)
        for (auto & e : entities_backup)
            e.model = this;

    // Pointing to the right place. This needs
    // to be done afterwards since the state zero is set as a function
    // of the population.
    db.model = this;
    db.user_data.model = this;

    if (use_queuing)
        queue.model = this;

    agents_data = model.agents_data;
    agents_data_ncols = model.agents_data_ncols;

    // Finally, seeds are resetted automatically based on the original
    // engine
    seed(
        static_cast<size_t>(
                std::floor(runif() * std::numeric_limits<size_t>::max())
            )
    );

}

template<typename TSeq>
inline Model<TSeq>::Model(Model<TSeq> && model) :
    name(std::move(model.name)),
    db(std::move(model.db)),
    population(std::move(model.population)),
    agents_data(std::move(model.agents_data)),
    agents_data_ncols(std::move(model.agents_data_ncols)),
    directed(std::move(model.directed)),
    // Virus
    viruses(std::move(model.viruses)),
    prevalence_virus(std::move(model.prevalence_virus)),
    prevalence_virus_as_proportion(std::move(model.prevalence_virus_as_proportion)),
    viruses_dist_funs(std::move(model.viruses_dist_funs)),
    // Tools
    tools(std::move(model.tools)),
    prevalence_tool(std::move(model.prevalence_tool)),
    prevalence_tool_as_proportion(std::move(model.prevalence_tool_as_proportion)),
    tools_dist_funs(std::move(model.tools_dist_funs)),
    // Entities
    entities(std::move(model.entities)),
    entities_backup(std::move(model.entities_backup)),
    // prevalence_entity(std::move(model.prevalence_entity)),
    // prevalence_entity_as_proportion(std::move(model.prevalence_entity_as_proportion)),
    // entities_dist_funs(std::move(model.entities_dist_funs)),
    // Pseudo-RNG
    engine(std::move(model.engine)),
    runifd(std::move(model.runifd)),
    rnormd(std::move(model.rnormd)),
    rgammad(std::move(model.rgammad)),
    rlognormald(std::move(model.rlognormald)),
    rexpd(std::move(model.rexpd)),
    // Rewiring
    rewire_fun(std::move(model.rewire_fun)),
    rewire_prop(std::move(model.rewire_prop)),
    parameters(std::move(model.parameters)),
    // Others
    ndays(model.ndays),
    pb(std::move(model.pb)),
    status_fun(std::move(model.status_fun)),
    states_labels(std::move(model.states_labels)),
    nstatus(model.nstatus),
    verbose(model.verbose),
    current_date(std::move(model.current_date)),
    global_action_functions(std::move(model.global_action_functions)),
    global_action_dates(std::move(model.global_action_dates)),
    queue(std::move(model.queue)),
    use_queuing(model.use_queuing),
    array_double_tmp(model.array_double_tmp.size()),
    array_virus_tmp(model.array_virus_tmp.size())
{

    db.model = this;
    db.user_data.model = this;

    if (use_queuing)
        queue.model = this;

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::operator=(const Model<TSeq> & m)
{
    name = m.name;

    population        = m.population;
    population_backup = m.population_backup;

    for (auto & p : population)
        p.model = this;

    if (population_backup.size() != 0)
        for (auto & p : population_backup)
            p.model = this;

    for (auto & e : entities)
        e.model = this;

    if (entities_backup.size() != 0)
        for (auto & e : entities_backup)
            e.model = this;

    db = m.db;
    db.model = this;
    db.user_data.model = this;

    directed = m.directed;
    
    viruses                        = m.viruses;
    prevalence_virus               = m.prevalence_virus;
    prevalence_virus_as_proportion = m.prevalence_virus_as_proportion;
    viruses_dist_funs              = m.viruses_dist_funs;

    tools                         = m.tools;
    prevalence_tool               = m.prevalence_tool;
    prevalence_tool_as_proportion = m.prevalence_tool_as_proportion;
    tools_dist_funs               = m.tools_dist_funs;
    
    entities        = m.entities;
    entities_backup = m.entities_backup;
    // prevalence_entity = m.prevalence_entity;
    // prevalence_entity_as_proportion = m.prevalence_entity_as_proportion;
    // entities_dist_funs = m.entities_dist_funs;
    
    rewire_fun  = m.rewire_fun;
    rewire_prop = m.rewire_prop;

    parameters = m.parameters;
    ndays      = m.ndays;
    pb         = m.pb;

    status_fun    = m.status_fun;
    states_labels = m.states_labels;
    nstatus       = m.nstatus;

    verbose     = m.verbose;

    current_date = m.current_date;

    global_action_functions = m.global_action_functions;
    global_action_dates     = m.global_action_dates;

    queue       = m.queue;
    use_queuing = m.use_queuing;

    // Making sure population is passed correctly
    // Pointing to the right place
    db.model = this;
    db.user_data.model = this;

    agents_data            = m.agents_data;
    agents_data_ncols = m.agents_data_ncols;

    // Figure out the queuing
    if (use_queuing)
        queue.model = this;

    // Finally, seeds are resetted automatically based on the original
    // engine
    seed(
        static_cast<size_t>(
                std::floor(runif() * std::numeric_limits<size_t>::max())
            )
    );

    array_double_tmp.resize(m.array_double_tmp.size());
    array_virus_tmp.resize(m.array_virus_tmp.size());

    return *this;

}

template<typename TSeq>
inline DataBase<TSeq> & Model<TSeq>::get_db()
{
    return db;
}

template<typename TSeq>
inline std::vector<Agent<TSeq>> & Model<TSeq>::get_agents()
{
    return population;
}

template<typename TSeq>
inline std::vector<Entity<TSeq>> & Model<TSeq>::get_entities()
{
    return entities;
}

template<typename TSeq>
inline void Model<TSeq>::agents_smallworld(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    bool d,
    epiworld_double p
)
{
    agents_from_adjlist(
        rgraph_smallworld(n, k, p, d, *this)
    );
}

template<typename TSeq>
inline void Model<TSeq>::agents_empty_graph(
    epiworld_fast_uint n
) 
{

    // Resizing the people
    population.clear();
    population.resize(n, Agent<TSeq>());

    // Filling the model and ids
    size_t i = 0u;
    for (auto & p : population)
    {
        p.id = i++;
        p.model = this;
    }
    

}

// template<typename TSeq>
// inline void Model<TSeq>::set_rand_engine(std::mt19937 & eng)
// {
//     engine = std::make_shared< std::mt19937 >(eng);
// }

template<typename TSeq>
inline void Model<TSeq>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::make_shared<std::gamma_distribution<>>(alpha,beta);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_norm(epiworld_double mean, epiworld_double sd)
{ 
    rnormd  = std::make_shared<std::normal_distribution<>>(mean, sd);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_unif(epiworld_double a, epiworld_double b)
{ 
    runifd  = std::make_shared<std::uniform_real_distribution<>>(a, b);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_lognormal(epiworld_double mean, epiworld_double shape)
{ 
    rlognormald  = std::make_shared<std::lognormal_distribution<>>(mean, shape);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_exp(epiworld_double lambda)
{ 
    rexpd  = std::make_shared<std::exponential_distribution<>>(lambda);
}

template<typename TSeq>
inline epiworld_double & Model<TSeq>::operator()(std::string pname) {

    if (parameters.find(pname) == parameters.end())
        throw std::range_error("The parameter "+ pname + "is not in the model.");

    return parameters[pname];

}

template<typename TSeq>
inline size_t Model<TSeq>::size() const {
    return population.size();
}

template<typename TSeq>
inline void Model<TSeq>::dist_virus()
{

    // Starting first infection
    int n = size();
    std::vector< size_t > idx(n);

    int n_left = n;
    std::iota(idx.begin(), idx.end(), 0);

    for (size_t v = 0u; v < viruses.size(); ++v)
    {

        if (viruses_dist_funs[v])
        {

            viruses_dist_funs[v](*viruses[v], this);

        } else {

            // Picking how many
            int nsampled;
            if (prevalence_virus_as_proportion[v])
            {
                nsampled = static_cast<int>(std::floor(prevalence_virus[v] * size()));
            }
            else
            {
                nsampled = static_cast<int>(prevalence_virus[v]);
            }

            if (nsampled > static_cast<int>(size()))
                throw std::range_error("There are only " + std::to_string(size()) + 
                " individuals in the population. Cannot add the virus to " + std::to_string(nsampled));


            VirusPtr<TSeq> virus = viruses[v];
            
            while (nsampled > 0)
            {

                int loc = static_cast<epiworld_fast_uint>(floor(runif() * (n_left--)));

                Agent<TSeq> & agent = population[idx[loc]];
                
                // Adding action
                agent.add_virus(
                    virus,
                    const_cast<Model<TSeq> * >(this),
                    virus->status_init,
                    virus->queue_init
                    );

                // Adjusting sample
                nsampled--;
                std::swap(idx[loc], idx[n_left]);

            }

        }

        // Apply the actions
        actions_run();
    }

}

template<typename TSeq>
inline void Model<TSeq>::dist_tools()
{

    // Starting first infection
    int n = size();
    std::vector< size_t > idx(n);
    for (epiworld_fast_uint t = 0; t < tools.size(); ++t)
    {

        if (tools_dist_funs[t])
        {

            tools_dist_funs[t](*tools[t], this);

        } else {

            // Picking how many
            int nsampled;
            if (prevalence_tool_as_proportion[t])
            {
                nsampled = static_cast<int>(std::floor(prevalence_tool[t] * size()));
            }
            else
            {
                nsampled = static_cast<int>(prevalence_tool[t]);
            }

            if (nsampled > static_cast<int>(size()))
                throw std::range_error("There are only " + std::to_string(size()) + 
                " individuals in the population. Cannot add the tool to " + std::to_string(nsampled));
            
            ToolPtr<TSeq> tool = tools[t];

            int n_left = n;
            std::iota(idx.begin(), idx.end(), 0);
            while (nsampled > 0)
            {
                int loc = static_cast<epiworld_fast_uint>(floor(runif() * n_left--));
                
                population[idx[loc]].add_tool(
                    tool,
                    const_cast< Model<TSeq> * >(this),
                    tool->status_init, tool->queue_init
                    );
                
                nsampled--;

                std::swap(idx[loc], idx[n_left]);

            }

        }

        // Apply the actions
        actions_run();

    }

}

// template<typename TSeq>
// inline void Model<TSeq>::dist_entities()
// {

//     // Starting first infection
//     int n = size();
//     std::vector< size_t > idx(n);
//     for (epiworld_fast_uint e = 0; e < entities.size(); ++e)
//     {

//         if (entities_dist_funs[e])
//         {

//             entities_dist_funs[e](entities[e], this);

//         } else {

//             // Picking how many
//             int nsampled;
//             if (prevalence_entity_as_proportion[e])
//             {
//                 nsampled = static_cast<int>(std::floor(prevalence_entity[e] * size()));
//             }
//             else
//             {
//                 nsampled = static_cast<int>(prevalence_entity[e]);
//             }

//             if (nsampled > static_cast<int>(size()))
//                 throw std::range_error("There are only " + std::to_string(size()) + 
//                 " individuals in the population. Cannot add the entity to " + std::to_string(nsampled));
            
//             Entity<TSeq> & entity = entities[e];

//             int n_left = n;
//             std::iota(idx.begin(), idx.end(), 0);
//             while (nsampled > 0)
//             {
//                 int loc = static_cast<epiworld_fast_uint>(floor(runif() * n_left--));
                
//                 population[idx[loc]].add_entity(entity, this, entity.status_init, entity.queue_init);
                
//                 nsampled--;

//                 std::swap(idx[loc], idx[n_left]);

//             }

//         }

//         // Apply the actions
//         actions_run();

//     }

// }

template<typename TSeq>
inline void Model<TSeq>::chrono_start() {
    time_start = std::chrono::steady_clock::now();
}

template<typename TSeq>
inline void Model<TSeq>::chrono_end() {
    time_end = std::chrono::steady_clock::now();
    time_elapsed += (time_end - time_start);
    n_replicates++;
}

template<typename TSeq>
inline void Model<TSeq>::set_backup()
{

    population_backup = population;

    entities_backup = entities;

}

// template<typename TSeq>
// inline void Model<TSeq>::restore_backup()
// {

//     // Restoring the data
//     population = *population_backup;
//     entities   = *entities_backup;

//     // And correcting the pointer
//     for (auto & p : population)
//         p.model = this;

//     for (auto & e : entities)
//         e.model = this;

// }

template<typename TSeq>
inline std::mt19937 & Model<TSeq>::get_rand_endgine()
{
    return engine;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif() {
    // CHECK_INIT()
    return runifd(engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif(epiworld_double a, epiworld_double b) {
    // CHECK_INIT()
    return runifd(engine) * (b - a) + a;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm() {
    // CHECK_INIT()
    return rnormd(engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm(epiworld_double mean, epiworld_double sd) {
    // CHECK_INIT()
    return rnormd(engine) * sd + mean;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma() {
    return rgammad(engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma(epiworld_double alpha, epiworld_double beta) {
    auto old_param = rgammad.param();
    rgammad.param(std::gamma_distribution<>::param_type(alpha, beta));
    epiworld_double ans = rgammad(engine);
    rgammad.param(old_param);
    return ans;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp() {
    return rexpd(engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp(epiworld_double lambda) {
    auto old_param = rexpd.param();
    rexpd.param(std::exponential_distribution<>::param_type(lambda));
    epiworld_double ans = rexpd(engine);
    rexpd.param(old_param);
    return ans;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal() {
    return rlognormald(engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal(epiworld_double mean, epiworld_double shape) {
    auto old_param = rlognormald.param();
    rlognormald.param(std::lognormal_distribution<>::param_type(mean, shape));
    epiworld_double ans = rlognormald(engine);
    rlognormald.param(old_param);
    return ans;
}

template<typename TSeq>
inline void Model<TSeq>::seed(size_t s) {
    this->engine.seed(s);
}

template<typename TSeq>
inline void Model<TSeq>::add_virus(Virus<TSeq> v, epiworld_double preval)
{

    if (preval > 1.0)
        throw std::range_error("Prevalence of virus cannot be above 1.0");

    if (preval < 0.0)
        throw std::range_error("Prevalence of virus cannot be negative");

    // Checking the state
    epiworld_fast_int init_, post_, rm_;
    v.get_state(&init_, &post_, &rm_);

    if (init_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -init- state."
            );
    else if (post_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -post- state."
            );
    
    // Recording the variant
    db.record_variant(v);

    // Adding new virus
    viruses.push_back(std::make_shared< Virus<TSeq> >(v));
    prevalence_virus.push_back(preval);
    prevalence_virus_as_proportion.push_back(true);
    viruses_dist_funs.push_back(nullptr);

}

template<typename TSeq>
inline void Model<TSeq>::add_virus_n(Virus<TSeq> v, epiworld_fast_uint preval)
{

    // Checking the ids
    epiworld_fast_int init_, post_, rm_;
    v.get_state(&init_, &post_, &rm_);

    if (init_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -init- state."
            );
    else if (post_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -post- state."
            );

    // Setting the id
    db.record_variant(v);

    // Adding new virus
    viruses.push_back(std::make_shared< Virus<TSeq> >(v));
    prevalence_virus.push_back(preval);
    prevalence_virus_as_proportion.push_back(false);
    viruses_dist_funs.push_back(nullptr);

}

template<typename TSeq>
inline void Model<TSeq>::add_virus_fun(Virus<TSeq> v, VirusToAgentFun<TSeq> fun)
{

    // Checking the ids
    epiworld_fast_int init_, post_, rm_;
    v.get_state(&init_, &post_, &rm_);

    if (init_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -init- state."
            );
    else if (post_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -post- state."
            );

    // Setting the id
    db.record_variant(v);
    // v.set_id(viruses.size());

    // Adding new virus
    viruses.push_back(std::make_shared< Virus<TSeq> >(v));
    prevalence_virus.push_back(0.0);
    prevalence_virus_as_proportion.push_back(false);
    viruses_dist_funs.push_back(fun);

}

template<typename TSeq>
inline void Model<TSeq>::add_tool(Tool<TSeq> t, epiworld_double preval)
{

    if (preval > 1.0)
        throw std::range_error("Prevalence of tool cannot be above 1.0");

    if (preval < 0.0)
        throw std::range_error("Prevalence of tool cannot be negative");

    db.record_tool(t);

    // Adding the tool to the model (and database.)
    tools.push_back(std::make_shared< Tool<TSeq> >(t));
    prevalence_tool.push_back(preval);
    prevalence_tool_as_proportion.push_back(true);
    tools_dist_funs.push_back(nullptr);

}

template<typename TSeq>
inline void Model<TSeq>::add_tool_n(Tool<TSeq> t, epiworld_fast_uint preval)
{
    
    db.record_tool(t);

    tools.push_back(std::make_shared<Tool<TSeq> >(t));
    prevalence_tool.push_back(preval);
    prevalence_tool_as_proportion.push_back(false);
    tools_dist_funs.push_back(nullptr);
}

template<typename TSeq>
inline void Model<TSeq>::add_tool_fun(Tool<TSeq> t, ToolToAgentFun<TSeq> fun)
{
    
    db.record_tool(t);
    
    tools.push_back(std::make_shared<Tool<TSeq> >(t));
    prevalence_tool.push_back(0.0);
    prevalence_tool_as_proportion.push_back(false);
    tools_dist_funs.push_back(fun);
}


template<typename TSeq>
inline void Model<TSeq>::add_entity(Entity<TSeq> e)
{

    e.model = this;
    e.id = entities.size();
    entities.push_back(e);

}

template<typename TSeq>
inline void Model<TSeq>::rm_virus(size_t virus_pos)
{

    if (viruses.size() <= virus_pos)
        throw std::range_error(
            std::string("The specified virus (") +
            std::to_string(virus_pos) +
            std::string(") is out of range. ") +
            std::string("There are only ") +
            std::to_string(viruses.size()) +
            std::string(" viruses.")
            );

    // Flipping with the last one
    std::swap(viruses[virus_pos], viruses[viruses.size() - 1]);
    std::swap(viruses_dist_funs[virus_pos], viruses_dist_funs[viruses.size() - 1]);
    std::swap(prevalence_virus[virus_pos], prevalence_virus[viruses.size() - 1]);
    std::swap(prevalence_virus_as_proportion[virus_pos], prevalence_virus_as_proportion[viruses.size() - 1]);

    viruses.pop_back();
    viruses_dist_funs.pop_back();
    prevalence_virus.pop_back();
    prevalence_virus_as_proportion.pop_back();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::rm_tool(size_t tool_pos)
{

    if (tools.size() <= tool_pos)
        throw std::range_error(
            std::string("The specified tool (") +
            std::to_string(tool_pos) +
            std::string(") is out of range. ") +
            std::string("There are only ") +
            std::to_string(tools.size()) +
            std::string(" tools.")
            );

    // Flipping with the last one
    std::swap(tools[tool_pos], tools[tools.size() - 1]);
    std::swap(tools_dist_funs[tool_pos], tools_dist_funs[tools.size() - 1]);
    std::swap(prevalence_tool[tool_pos], prevalence_tool[tools.size() - 1]);
    std::swap(prevalence_tool_as_proportion[tool_pos], prevalence_tool_as_proportion[tools.size() - 1]);


    tools.pop_back();
    tools_dist_funs.pop_back();
    prevalence_tool.pop_back();
    prevalence_tool_as_proportion.pop_back();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::load_agents_entities_ties(
    std::string fn,
    int skip
    )
{

    int i,j;
    std::ifstream filei(fn);

    if (!filei)
        throw std::logic_error("The file " + fn + " was not found.");

    int linenum = 0;
    std::vector< epiworld_fast_uint > source_;
    std::vector< std::vector< epiworld_fast_uint > > target_(entities.size());

    target_.reserve(1e5);

    while (!filei.eof())
    {

        if (linenum++ < skip)
            continue;

        filei >> i >> j;

        // Looking for exceptions
        if (filei.bad())
            throw std::logic_error(
                "I/O error while reading the file " +
                fn
            );

        if (filei.fail())
            break;

        if (i >= static_cast<int>(this->size()))
            throw std::range_error(
                "The agent["+std::to_string(linenum)+"] = " + std::to_string(i) +
                " is above the max id " + std::to_string(this->size() - 1)
                );

        if (j >= static_cast<int>(this->entities.size()))
            throw std::range_error(
                "The entity["+std::to_string(linenum)+"] = " + std::to_string(j) +
                " is above the max id " + std::to_string(this->entities.size() - 1)
                );

        target_[j].push_back(i);

        population[i].add_entity(entities[j], nullptr);

    }

    // // Iterating over entities
    // for (size_t e = 0u; e < entities.size(); ++e)
    // {

    //     // This entity will have individuals assigned to it, so we add it
    //     if (target_[e].size() > 0u)
    //     {

    //         // Filling in the gaps
    //         prevalence_entity[e] = static_cast<epiworld_double>(target_[e].size());
    //         prevalence_entity_as_proportion[e] = false;

    //         // Generating the assignment function
    //         auto who = target_[e];
    //         entities_dist_funs[e] =
    //             [who](Entity<TSeq> & e, Model<TSeq>* m) -> void {

    //                 for (auto w : who)
    //                     m->population[w].add_entity(e, m, e.status_init, e.queue_init);
                    
    //                 return;
                    
    //             };

    //     }

    // }

    return;

}

template<typename TSeq>
inline void Model<TSeq>::agents_from_adjlist(
    std::string fn,
    int size,
    int skip,
    bool directed
    ) {

    AdjList al;
    al.read_edgelist(fn, size, skip, directed);
    this->agents_from_adjlist(al);

}

template<typename TSeq>
inline void Model<TSeq>::agents_from_edgelist(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) {

    AdjList al(source, target, size, directed);
    agents_from_adjlist(al);

}

template<typename TSeq>
inline void Model<TSeq>::agents_from_adjlist(AdjList al) {

    // Resizing the people
    agents_empty_graph(al.vcount());
    
    const auto & tmpdat = al.get_dat();
    
    for (size_t i = 0u; i < tmpdat.size(); ++i)
    {

        // population[i].id    = i;
        population[i].model = this;

        for (const auto & link: tmpdat[i])
        {

            population[i].add_neighbor(
                population[link.first],
                true, true
                );

        }

    }

    #ifdef EPI_DEBUG
    for (auto & p: population)
    {
        if (p.id >= static_cast<int>(al.vcount()))
            throw std::logic_error(
                "Agent's id cannot be negative above or equal to the number of agents!");
    }
    #endif

}

template<typename TSeq>
inline bool Model<TSeq>::is_directed() const
{
    if (population.size() == 0u)
        throw std::logic_error("The population hasn't been initialized.");

    return directed;
}

template<typename TSeq>
inline int Model<TSeq>::today() const {

    if (ndays == 0)
      return 0;

    return this->current_date;
}

template<typename TSeq>
inline void Model<TSeq>::next() {

    db.record();
    ++this->current_date;
    
    // Advancing the progress bar
    if ((this->current_date >= 1) && verbose)
        pb.next();

    return ;
}

template<typename TSeq>
inline void Model<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
) 
{

    if (size() == 0u)
        throw std::logic_error("There's no agents in this model!");

    if (nstatus == 0u)
        throw std::logic_error(
            std::string("No states registered in this model. ") +
            std::string("At least one state should be included. See the function -Model::add_state()-")
            );

    // Setting up the number of steps
    this->ndays = ndays;

    if (seed >= 0)
        engine.seed(seed);

    array_double_tmp.resize(size()/2, 0.0);
    array_virus_tmp.resize(size()/2);

    // Checking whether the proposed state in/out/removed
    // are valid
    epiworld_fast_int _init, _end, _removed;
    int nstatus_int = static_cast<int>(nstatus);
    for (auto & v : viruses)
    {
        v->get_state(&_init, &_end, &_removed);
        
        // Negative unspecified state
        if (((_init != -99) && (_init < 0)) || (_init >= nstatus_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstatus - 1));

        // Negative unspecified state
        if (((_end != -99) && (_end < 0)) || (_end >= nstatus_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstatus - 1));

        if (((_removed != -99) && (_removed < 0)) || (_removed >= nstatus_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstatus - 1));

    }

    for (auto & t : tools)
    {
        t->get_state(&_init, &_end);
        
        // Negative unspecified state
        if (((_init != -99) && (_init < 0)) || (_init >= nstatus_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstatus - 1));

        // Negative unspecified state
        if (((_end != -99) && (_end < 0)) || (_end >= nstatus_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstatus - 1));

    }

    // Starting first infection and tools
    reset();

    // Initializing the simulation
    chrono_start();
    EPIWORLD_RUN((*this))
    {

        #ifdef EPI_DEBUG
        db.n_transmissions_potential = 0;
        db.n_transmissions_today = 0;
        #endif

        // We can execute these components in whatever order the
        // user needs.
        this->update_state();
    
        // We start with the global actions
        this->run_global_actions();

        // In this case we are applying degree sequence rewiring
        // to change the network just a bit.
        this->rewire();

        // This locks all the changes
        this->next();

        // Mutation must happen at the very end of all
        this->mutate_variant();

    }

    // The last reaches the end...
    this->current_date--;

    chrono_end();

}

template<typename TSeq>
inline void Model<TSeq>::run_multiple(
    epiworld_fast_uint ndays,
    epiworld_fast_uint nexperiments,
    int seed_,
    std::function<void(size_t,Model<TSeq>*)> fun,
    bool reset,
    bool verbose,
    int nthreads
)
{

    if (seed_ >= 0)
        this->seed(seed_);

    // Seeds will be reproducible by default
    std::vector< int > seeds_n(nexperiments);
    #ifdef EPI_DEBUG
    std::fill(
        seeds_n.begin(),
        seeds_n.end(),
        std::floor(
            runif() * static_cast<double>(std::numeric_limits<int>::max())
        )
        );
    #else
    for (auto & s : seeds_n)
    {
        s = static_cast<int>(
            std::floor(
                runif() * static_cast<double>(std::numeric_limits<int>::max())
                )
        );
    }
    #endif

    EPI_DEBUG_NOTIFY_ACTIVE()

    bool old_verb = this->verbose;
    verbose_off();

    // Setting up backup
    if (reset)
        set_backup();

    #ifdef _OPENMP

    omp_set_num_threads(nthreads);

    // Generating copies of the model
    std::vector< Model<TSeq> * > these;
    for (size_t i = 0; i < static_cast<size_t>(std::max(nthreads - 1, 0)); ++i)
        these.push_back(clone_ptr());

    // Figuring out how many replicates
    std::vector< size_t > nreplicates(nthreads, 0);
    std::vector< size_t > nreplicates_csum(nthreads, 0);
    size_t sums = 0u;
    for (int i = 0; i < nthreads; ++i)
    {
        nreplicates[i] = static_cast<epiworld_fast_uint>(
            std::floor(nexperiments/nthreads)
            );
        
        // This takes the cumsum
        nreplicates_csum[i] = sums;

        sums += nreplicates[i];

    }

    if (sums < nexperiments)
        nreplicates[nthreads - 1] += (nexperiments - sums);

    Progress pb_multiple(
        nreplicates[0u],
        EPIWORLD_PROGRESS_BAR_WIDTH
        );

    if (verbose)
    {

        printf_epiworld(
            "Starting multiple runs (%i) using %i thread(s)\n", 
            static_cast<int>(nexperiments),
            static_cast<int>(nthreads)
        );

        pb_multiple.start();

    }

    #pragma omp parallel shared(these, nreplicates, nreplicates_csum, seeds_n) \
        firstprivate(nexperiments, nthreads, fun, reset, verbose, pb_multiple, ndays) \
        default(shared)
    {

        auto iam = omp_get_thread_num();

        for (size_t n = 0u; n < nreplicates[iam]; ++n)
        {
            
            if (iam == 0)
            {

                // Initializing the seed
                run(ndays, seeds_n[n]);

                if (fun)
                    fun(n, this);

                // Only the first one prints
                if (verbose)
                    pb_multiple.next();

            } else {

                // Initializing the seed
                these[iam - 1]->run(ndays, seeds_n[nreplicates_csum[iam] + n]);

                if (fun)
                    fun(
                        n + nreplicates_csum[iam],
                        these[iam - 1]
                        );

            }

        }
        
    }

    // Adjusting the number of replicates
    n_replicates += (nexperiments - nreplicates[0u]);

    for (auto & ptr : these)
        delete ptr;

    #else
    if (reset)
        set_backup();

    Progress pb_multiple(
        nexperiments,
        EPIWORLD_PROGRESS_BAR_WIDTH
        )
        ;
    if (verbose)
    {

        printf_epiworld(
            "Starting multiple runs (%i)\n", 
            static_cast<int>(nexperiments)
        );

        pb_multiple.start();

    }

    for (size_t n = 0u; n < nexperiments; ++n)
    {

        run(ndays, seeds_n[n]);

        if (fun)
            fun(n, this);

        if (verbose)
            pb_multiple.next();
    
    }
    #endif

    if (verbose)
        pb_multiple.end();

    if (old_verb)
        verbose_on();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::update_state() {

    // Next state
    if (use_queuing)
    {
        int i = -1;
        for (auto & p: population)
            if (queue[++i] > 0)
            {
                if (status_fun[p.state])
                    status_fun[p.state](&p, this);
            }

    }
    else
    {

        for (auto & p: population)
            if (status_fun[p.state])
                    status_fun[p.state](&p, this);

    }

    actions_run();
    
}



template<typename TSeq>
inline void Model<TSeq>::mutate_variant() {

    if (use_queuing)
    {

        int i = -1;
        for (auto & p: population)
        {

            if (queue[++i] == 0)
                continue;

            if (p.n_viruses > 0u)
                for (auto & v : p.get_viruses())
                    v->mutate(this);

        }

    }
    else 
    {

        for (auto & p: population)
        {

            if (p.n_viruses > 0u)
                for (auto & v : p.get_viruses())
                    v->mutate(this);

        }

    }
    

}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_variants() const {
    return db.size();
}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_tools() const {
    return tools.size();
}

template<typename TSeq>
inline epiworld_fast_uint Model<TSeq>::get_ndays() const {
    return ndays;
}

template<typename TSeq>
inline epiworld_fast_uint Model<TSeq>::get_n_replicates() const
{
    return n_replicates;
}

template<typename TSeq>
inline void Model<TSeq>::set_ndays(epiworld_fast_uint ndays) {
    this->ndays = ndays;
}

template<typename TSeq>
inline bool Model<TSeq>::get_verbose() const {
    return verbose;
}

template<typename TSeq>
inline void Model<TSeq>::verbose_on() {
    verbose = true;
}

template<typename TSeq>
inline void Model<TSeq>::verbose_off() {
    verbose = false;
}

template<typename TSeq>
inline void Model<TSeq>::set_rewire_fun(
    std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> fun
    ) {
    rewire_fun = fun;
}

template<typename TSeq>
inline void Model<TSeq>::set_rewire_prop(epiworld_double prop)
{

    if (prop < 0.0)
        throw std::range_error("Proportions cannot be negative.");

    if (prop > 1.0)
        throw std::range_error("Proportions cannot be above 1.0.");

    rewire_prop = prop;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::get_rewire_prop() const {
    return rewire_prop;
}

template<typename TSeq>
inline void Model<TSeq>::rewire() {

    if (rewire_fun)
        rewire_fun(&population, this, rewire_prop);
}


template<typename TSeq>
inline void Model<TSeq>::write_data(
    std::string fn_variant_info,
    std::string fn_variant_hist,
    std::string fn_tool_info,
    std::string fn_tool_hist,
    std::string fn_total_hist,
    std::string fn_transmission,
    std::string fn_transition,
    std::string fn_reproductive_number,
    std::string fn_generation_time
    ) const
{

    db.write_data(
        fn_variant_info, fn_variant_hist,
        fn_tool_info, fn_tool_hist,
        fn_total_hist, fn_transmission, fn_transition,
        fn_reproductive_number, fn_generation_time
        );

}

template<typename TSeq>
inline void Model<TSeq>::write_edgelist(
    std::string fn
    ) const
{

    // Figuring out the writing sequence
    std::vector< const Agent<TSeq> * > wseq(size());
    for (const auto & p: population)
        wseq[p.id] = &p;

    std::ofstream efile(fn, std::ios_base::out);
    efile << "source target\n";
    if (this->is_directed())
    {

        for (const auto & p : wseq)
        {
            for (auto & n : p->neighbors)
                efile << p->id << " " << n << "\n";
        }

    } else {

        for (const auto & p : wseq)
        {
            for (auto & n : p->neighbors)
                if (static_cast<int>(p->id) <= static_cast<int>(n))
                    efile << p->id << " " << n << "\n";
        }

    }

}

template<typename TSeq>
inline std::map<std::string,epiworld_double> & Model<TSeq>::params()
{
    return parameters;
}

template<typename TSeq>
inline void Model<TSeq>::reset() {

    // Restablishing people
    pb = Progress(ndays, 80);

    if (population_backup.size() != 0u)
    {
        population = population_backup;

        #ifdef EPI_DEBUG
        for (size_t i = 0; i < population.size(); ++i)
        {

            if (population[i] != population_backup[i])
                throw std::logic_error("Model::reset population doesn't match.");

        }
        #endif

    } 
    else
    {
        for (auto & p : population)
            p.reset();
    }
        
    if (entities_backup.size() != 0)
    {
        entities = entities_backup;

        #ifdef EPI_DEBUG
        for (size_t i = 0; i < entities.size(); ++i)
        {

            if (entities[i] != entities_backup[i])
                throw std::logic_error("Model::reset entities don't match.");

        }
        #endif
        
    }
    else
    {
        for (auto & e: entities)
            e.reset();
    }
    
    current_date = 0;

    db.reset();

    // This also clears the queue
    if (use_queuing)
        queue.reset();

    // Re distributing tools and virus
    dist_virus();
    dist_tools();

    // Recording the original state (at time 0) and advancing
    // to time 1
    next();


}

// Too big to keep here
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//model-meat-print.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODEL_MEAT_PRINT_HPP
#define EPIWORLD_MODEL_MEAT_PRINT_HPP

template<typename TSeq>
inline void Model<TSeq>::print(bool lite) const
{

    if (lite)
    {
        size_t nchar = 0u;
        std::string fmt = " - %-" + std::to_string(nchar + 1) + "s: ";

        for (auto & p : states_labels)
            if (p.length() > nchar)
                nchar = p.length();

        if (today() != 0)
            fmt = "  - (%" + std::to_string(nstatus).length() +
                std::string("d) %-") + std::to_string(nchar) + "s : %" +
                std::to_string(std::to_string(size()).length()) + "i -> %i\n";
        else
            fmt = "  - (%" + std::to_string(nstatus).length() +
                std::string("d) %-") + std::to_string(nchar) + "s : %i\n";

        printf_epiworld("\nDistribution of the population at time %i:\n", today());
        for (size_t s = 0u; s < nstatus; ++s)
        {
            if (today() != 0)
            {

                printf_epiworld(
                    fmt.c_str(),
                    s,
                    states_labels[s].c_str(),
                    db.hist_total_counts[s],
                    db.today_total[ s ]
                    );

            }
            else
            {

                printf_epiworld(
                    fmt.c_str(),
                    s,
                    states_labels[s].c_str(),
                    db.today_total[ s ]
                    );

            }
        }

        return;
    }

    // Horizontal line
    std::string line = "";
    for (epiworld_fast_uint i = 0u; i < 80u; ++i)
        line += "_";

    // Prints a message if debugging is on
    EPI_DEBUG_NOTIFY_ACTIVE()

    printf_epiworld("\n%s\n%s\n\n",line.c_str(), "SIMULATION STUDY");

    printf_epiworld("Name of the model   : %s\n", (this->name == "") ? std::string("(none)").c_str() : name.c_str());
    printf_epiworld("Population size     : %i\n", static_cast<int>(size()));
    printf_epiworld("Number of entities  : %i\n", static_cast<int>(entities.size()));
    printf_epiworld("Days (duration)     : %i (of %i)\n", today(), static_cast<int>(ndays));
    printf_epiworld("Number of variants  : %i\n", static_cast<int>(db.get_n_variants()));
    if (n_replicates > 0u)
    {
        std::string abbr;
        epiworld_double elapsed;
        epiworld_double total;
        get_elapsed("auto", &elapsed, &total, &abbr, false);
        printf_epiworld("Last run elapsed t  : %.2f%s\n", elapsed, abbr.c_str());
        if (n_replicates > 1u)
        {
            printf_epiworld("Total elapsed t     : %.2f%s (%i runs)\n", total, abbr.c_str(), static_cast<int>(n_replicates));
        }

        // Elapsed time in speed
        get_elapsed("microseconds", &elapsed, &total, &abbr, false);
        printf_epiworld("Last run speed      : %.2f million agents x day / second\n",
            static_cast<double>(this->size()) *
            static_cast<double>(this->get_ndays()) /
            static_cast<double>(elapsed)
            );
        if (n_replicates > 1u)
        {
            printf_epiworld("Average run speed   : %.2f million agents x day / second\n",
                static_cast<double>(this->size()) *
                static_cast<double>(this->get_ndays()) *
                static_cast<double>(n_replicates) /
                static_cast<double>(total)
            );
        }

    } else {
        printf_epiworld("Last run elapsed t  : -\n");
    }
    
    
    if (rewire_fun)
    {
        printf_epiworld("Rewiring            : on (%.2f)\n\n", rewire_prop);
    } else {
        printf_epiworld("Rewiring            : off\n\n");
    }
    

    printf_epiworld("Virus(es):\n");
    size_t n_variants_model = viruses.size();
    for (size_t i = 0u; i < n_variants_model; ++i)
    {    

        if ((n_variants_model > 10) && (i >= 10))
        {
            printf_epiworld(" ...and %li more variants...\n", n_variants_model - i);
            break;
        }

        if (i < n_variants_model)
        {

            if (prevalence_virus_as_proportion[i])
            {

                printf_epiworld(
                    " - %s (baseline prevalence: %.2f%%)\n",
                    viruses[i]->get_name().c_str(),
                    prevalence_virus[i] * 100.00
                );

            }
            else
            {

                printf_epiworld(
                    " - %s (baseline prevalence: %i seeds)\n",
                    viruses[i]->get_name().c_str(),
                    static_cast<int>(prevalence_virus[i])
                );

            }

        } else {

            printf_epiworld(
                " - %s (originated in the model...)\n",
                viruses[i]->get_name().c_str()
            );

        }

    }

    if (viruses.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    printf_epiworld("\nTool(s):\n");
    size_t n_tools_model = tools.size();
    for (size_t i = 0u; i < tools.size(); ++i)
    {   

        if ((n_tools_model > 10) && (i >= 10))
        {
            printf_epiworld(" ...and %li more tools...\n", n_tools_model - i);
            break;
        }

        if (i < n_tools_model)
        {
            if (prevalence_tool_as_proportion[i])
            {

                printf_epiworld(
                    " - %s (baseline prevalence: %.2f%%)\n",
                    tools[i]->get_name().c_str(),
                    prevalence_tool[i] * 100.0
                    );

            }
            else
            {

                printf_epiworld(
                    " - %s (baseline prevalence: %i seeds)\n",
                    tools[i]->get_name().c_str(),
                    static_cast<int>(prevalence_tool[i])
                    );

            }

        } else {

            printf_epiworld(
                " - %s (originated in the model...)\n",
                tools[i]->get_name().c_str()
            );

        }
        

    }

    if (tools.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    // Information about the parameters included
    printf_epiworld("\nModel parameters:\n");
    epiworld_fast_uint nchar = 0u;
    for (auto & p : parameters)
        if (p.first.length() > nchar)
            nchar = p.first.length();

    std::string fmt = " - %-" + std::to_string(nchar + 1) + "s: ";
    for (auto & p : parameters)
    {
        std::string fmt_tmp = fmt;
        if (std::fabs(p.second) < 0.0001)
            fmt_tmp += "%.1e\n";
        else
            fmt_tmp += "%.4f\n";

        printf_epiworld(
            fmt_tmp.c_str(),
            p.first.c_str(),
            p.second
        );
        
    }

    if (parameters.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    nchar = 0u;
    for (auto & p : states_labels)
        if (p.length() > nchar)
            nchar = p.length();

    

    if (today() != 0)
        fmt = "  - (%" + std::to_string(nstatus).length() +
            std::string("d) %-") + std::to_string(nchar) + "s : %" +
            std::to_string(std::to_string(size()).length()) + "i -> %i\n";
    else
        fmt = "  - (%" + std::to_string(nstatus).length() +
            std::string("d) %-") + std::to_string(nchar) + "s : %i\n";
        
    if (today() != 0)
    {
        printf_epiworld("\nDistribution of the population at time %i:\n", today());
        for (size_t s = 0u; s < nstatus; ++s)
        {

                printf_epiworld(
                    fmt.c_str(),
                    s,
                    states_labels[s].c_str(),
                    db.hist_total_counts[s],
                    db.today_total[ s ]
                    );

        }
            // else
            // {

            //     printf_epiworld(
            //         fmt.c_str(),
            //         s,
            //         states_labels[s].c_str(),
            //         db.today_total[ s ]
            //         );

            // }
    }

    if (today() != 0)
        (void) db.transition_probability(true);

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//model-meat-print.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



template<typename TSeq>
inline Model<TSeq> && Model<TSeq>::clone() const {

    // Step 1: Regen the individuals and make sure that:
    //  - Neighbors point to the right place
    //  - DB is pointing to the right place
    Model<TSeq> res(*this);

    // Removing old neighbors
    for (auto & p: res.population)
        p.neighbors.clear();
    
    // Rechecking individuals
    for (epiworld_fast_uint p = 0u; p < size(); ++p)
    {
        // Making room
        const Agent<TSeq> & agent_this = population[p];
        Agent<TSeq> & agent_res  = res.population[p];

        // Agent pointing to the right model and agent
        agent_res.model         = &res;
        agent_res.viruses.agent = &agent_res;
        agent_res.tools.agent   = &agent_res;

        // Readding
        std::vector< Agent<TSeq> * > neigh = agent_this.neighbors;
        for (epiworld_fast_uint n = 0u; n < neigh.size(); ++n)
        {
            // Point to the right neighbors
            int loc = res.population_ids[neigh[n]->get_id()];
            agent_res.add_neighbor(res.population[loc], true, true);

        }

    }

    return res;

}



template<typename TSeq>
inline void Model<TSeq>::add_state(
    std::string lab, 
    UpdateFun<TSeq> fun
)
{

    // Checking it doesn't match
    for (auto & s : states_labels)
        if (s == lab)
            throw std::logic_error("state \"" + s + "\" already registered.");

    states_labels.push_back(lab);
    status_fun.push_back(fun);
    nstatus++;

}


template<typename TSeq>
inline const std::vector< std::string > &
Model<TSeq>::get_state() const
{
    return states_labels;
}

template<typename TSeq>
inline const std::vector< UpdateFun<TSeq> > &
Model<TSeq>::get_state_fun() const
{
    return status_fun;
}

template<typename TSeq>
inline void Model<TSeq>::print_state_codes() const
{

    // Horizontal line
    std::string line = "";
    for (epiworld_fast_uint i = 0u; i < 80u; ++i)
        line += "_";

    printf_epiworld("\n%s\nSTATUS CODES\n\n", line.c_str());

    epiworld_fast_uint nchar = 0u;
    for (auto & p : states_labels)
        if (p.length() > nchar)
            nchar = p.length();
    
    std::string fmt = " %2i = %-" + std::to_string(nchar + 1 + 4) + "s\n";
    for (epiworld_fast_uint i = 0u; i < nstatus; ++i)
    {

        printf_epiworld(
            fmt.c_str(),
            i,
            (states_labels[i] + " (S)").c_str()
        );

    }

}



template<typename TSeq>
inline epiworld_double Model<TSeq>::add_param(
    epiworld_double initial_value,
    std::string pname
    ) {

    if (parameters.find(pname) == parameters.end())
        parameters[pname] = initial_value;
    
    return initial_value;

}

template<typename TSeq>
inline void Model<TSeq>::read_params(std::string fn)
{

    std::ifstream paramsfile(fn);

    if (!paramsfile)
        throw std::logic_error("The file " + fn + " was not found.");

    std::regex pattern("^([^:]+)\\s*[:]\\s*([0-9]+|[0-9]*\\.[0-9]+)?\\s*$");

    std::string line;
    std::smatch match;
    auto empty = std::sregex_iterator();

    while (std::getline(paramsfile, line))
    {

        // Is it a comment or an empty line?
        if (std::regex_match(line, std::regex("^([*].+|//.+|#.+|\\s*)$")))
            continue;

        // Finding the patter, if it doesn't match, then error
        std::regex_match(line, match, pattern);

        if (match.empty())
            throw std::logic_error("The line does not match parameters:\n" + line);

        // Capturing the number
        std::string anumber = match[2u].str() + match[3u].str();
        epiworld_double tmp_num = static_cast<epiworld_double>(
            std::strtod(anumber.c_str(), nullptr)
            );

        add_param(
            tmp_num,
            std::regex_replace(
                match[1u].str(),
                std::regex("^\\s+|\\s+$"),
                "")
        );

    }

}

template<typename TSeq>
inline epiworld_double Model<TSeq>::get_param(std::string pname)
{
    if (parameters.find(pname) == parameters.end())
        throw std::logic_error("The parameter " + pname + " does not exists.");

    return parameters[pname];
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::par(std::string pname)
{
    return parameters[pname];
}

#define DURCAST(tunit,txtunit) {\
        elapsed       = std::chrono::duration_cast<std::chrono:: tunit>(\
            time_end - time_start).count(); \
        elapsed_total = std::chrono::duration_cast<std::chrono:: tunit>(time_elapsed).count(); \
        abbr_unit     = txtunit;}

template<typename TSeq>
inline void Model<TSeq>::get_elapsed(
    std::string unit,
    epiworld_double * last_elapsed,
    epiworld_double * total_elapsed,
    std::string * unit_abbr,
    bool print
) const {

    // Preparing the result
    epiworld_double elapsed, elapsed_total;
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
    if (total_elapsed != nullptr)
        *total_elapsed = elapsed_total;
    if (unit_abbr != nullptr)
        *unit_abbr = abbr_unit;

    if (!print)
        return;

    if (n_replicates > 1u)
    {
        printf_epiworld("last run elapsed time : %.2f%s\n",
            elapsed, abbr_unit.c_str());
        printf_epiworld("total elapsed time    : %.2f%s\n",
            elapsed_total, abbr_unit.c_str());
        printf_epiworld("total runs            : %i\n",
            static_cast<int>(n_replicates));
        printf_epiworld("mean run elapsed time : %.2f%s\n",
            elapsed_total/static_cast<epiworld_double>(n_replicates), abbr_unit.c_str());

    } else {
        printf_epiworld("last run elapsed time : %.2f%s.\n", elapsed, abbr_unit.c_str());
    }
}

template<typename TSeq>
inline void Model<TSeq>::set_user_data(std::vector< std::string > names)
{
    db.set_user_data(names);
}

template<typename TSeq>
inline void Model<TSeq>::add_user_data(epiworld_fast_uint j, epiworld_double x)
{
    db.add_user_data(j, x);
}

template<typename TSeq>
inline void Model<TSeq>::add_user_data(std::vector<epiworld_double> x)
{
    db.add_user_data(x);
}

template<typename TSeq>
inline UserData<TSeq> & Model<TSeq>::get_user_data()
{
    return db.get_user_data();
}

template<typename TSeq>
inline void Model<TSeq>::add_global_action(
    std::function<void(Model<TSeq>*)> fun,
    int date
)
{

    global_action_functions.push_back(fun);
    global_action_dates.push_back(date);

}

template<typename TSeq>
inline void Model<TSeq>::run_global_actions()
{

    for (epiworld_fast_uint i = 0u; i < global_action_dates.size(); ++i)
    {

        if (global_action_dates[i] < 0)
        {

            global_action_functions[i](this);

        }
        else if (global_action_dates[i] == today())
        {

            global_action_functions[i](this);

        }

        actions_run();

    }

}

template<typename TSeq>
inline void Model<TSeq>::queuing_on()
{
    use_queuing = true;
}

template<typename TSeq>
inline void Model<TSeq>::queuing_off()
{
    use_queuing = false;
}

template<typename TSeq>
inline bool Model<TSeq>::is_queuing_on() const
{
    return use_queuing;
}

template<typename TSeq>
inline Queue<TSeq> & Model<TSeq>::get_queue()
{
    return queue;
}

template<typename TSeq>
inline const std::vector< VirusPtr<TSeq> > & Model<TSeq>::get_viruses() const
{
    return viruses;
}

template<typename TSeq>
const std::vector< ToolPtr<TSeq> > & Model<TSeq>::get_tools() const
{
    return tools;
}

template<typename TSeq>
inline void Model<TSeq>::set_agents_data(double * data_, size_t ncols_)
{
    agents_data = data_;
    agents_data_ncols = ncols_;
}

template<typename TSeq>
inline double * Model<TSeq>::get_agents_data() {
    return this->agents_data;
}

template<typename TSeq>
inline size_t Model<TSeq>::get_agents_data_ncols()  {
    return this->agents_data_ncols;
}


template<typename TSeq>
inline void Model<TSeq>::set_name(std::string name)
{
    this->name = name;
}

template<typename TSeq>
inline std::string Model<TSeq>::get_name() const 
{
    return this->name;
}

#define VECT_MATCH(a, b, c) \
    EPI_DEBUG_FAIL_AT_TRUE(a.size() != b.size(), c) \
    for (size_t __i = 0u; __i < a.size(); ++__i) \
    {\
        EPI_DEBUG_FAIL_AT_TRUE(a[__i] != b[__i], c) \
    }

template<typename TSeq>
inline bool Model<TSeq>::operator==(const Model<TSeq> & other) const
{
    EPI_DEBUG_FAIL_AT_TRUE(name != other.name, "names don't match")
    EPI_DEBUG_FAIL_AT_TRUE(db != other.db, "database don't match")

    VECT_MATCH(population, other.population, "population doesn't match")

    EPI_DEBUG_FAIL_AT_TRUE(
        using_backup != other.using_backup,
        "Model:: using_backup don't match"
        )
    
    if ((population_backup.size() != 0) & (other.population_backup.size() != 0))
    {
        for (size_t i = 0u; i < population_backup.size(); ++i)
        {
            if (population_backup[i] != other.population_backup[i])
                return false;
        }
        
    } else if ((population_backup.size() == 0) & (other.population_backup.size() != 0)) {
        return false;
    } else if ((population_backup.size() != 0) & (other.population_backup.size() == 0))
    {
        return false;
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        agents_data != other.agents_data,
        "Model:: agents_data don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        agents_data_ncols != other.agents_data_ncols,
        "Model:: agents_data_ncols don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        directed != other.directed,
        "Model:: directed don't match"
    )
    
    // Viruses -----------------------------------------------------------------
    EPI_DEBUG_FAIL_AT_TRUE(
        viruses.size() != other.viruses.size(),
        "Model:: viruses.size() don't match"
        )

    for (size_t i = 0u; i < viruses.size(); ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *viruses[i] != *other.viruses[i],
            "Model:: *viruses[i] don't match"
        )
            
    }

    VECT_MATCH(
        prevalence_virus,
        other.prevalence_virus,
        "virus prevalence don't match"
    )

    VECT_MATCH(
        prevalence_virus_as_proportion,
        other.prevalence_virus_as_proportion,
        "virus prevalence as prop don't match"
    )
    
    // Tools -------------------------------------------------------------------
    EPI_DEBUG_FAIL_AT_TRUE(
        tools.size() != other.tools.size(),
        "Model:: tools.size() don't match"
        )
        
    for (size_t i = 0u; i < tools.size(); ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *tools[i] != *other.tools[i],
            "Model:: *tools[i] don't match"
        )
            
    }

    VECT_MATCH(
        prevalence_tool, 
        other.prevalence_tool, 
        "tools prevalence don't match"
    )

    VECT_MATCH(
        prevalence_tool_as_proportion, 
        other.prevalence_tool_as_proportion, 
        "tools as prop don't match"
    )
    
    VECT_MATCH(
        entities,
        other.entities,
        "entities don't match"
    )

    if ((entities_backup.size() != 0) & (other.entities_backup.size() != 0))
    {
        
        for (size_t i = 0u; i < entities_backup.size(); ++i)
        {

            EPI_DEBUG_FAIL_AT_TRUE(
                entities_backup[i] != other.entities_backup[i],
                "Model:: entities_backup[i] don't match"
            )

        }
        
    } else if ((entities_backup.size() == 0) & (other.entities_backup.size() != 0)) {
        EPI_DEBUG_FAIL_AT_TRUE(true, "entities_backup don't match")
    } else if ((entities_backup.size() != 0) & (other.entities_backup.size() == 0))
    {
        EPI_DEBUG_FAIL_AT_TRUE(true, "entities_backup don't match")
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        rewire_prop != other.rewire_prop,
        "Model:: rewire_prop don't match"
    )
        
    EPI_DEBUG_FAIL_AT_TRUE(
        parameters.size() != other.parameters.size(),
        "Model:: () don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        parameters != other.parameters,
        "Model:: parameters don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        ndays != other.ndays,
        "Model:: ndays don't match"
    )
    
    VECT_MATCH(
        states_labels,
        other.states_labels,
        "state labels don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        nstatus != other.nstatus,
        "Model:: nstatus don't match"
    )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        verbose != other.verbose,
        "Model:: verbose don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        current_date != other.current_date,
        "Model:: current_date don't match"
    )

    VECT_MATCH(global_action_dates, other.global_action_dates, "global action dates don't match");

    EPI_DEBUG_FAIL_AT_TRUE(
        queue != other.queue,
        "Model:: queue don't match"
    )
    

    EPI_DEBUG_FAIL_AT_TRUE(
        use_queuing != other.use_queuing,
        "Model:: use_queuing don't match"
    )
    
    return true;

}

#undef VECT_MATCH
#undef DURCAST
#undef CASES_PAR
#undef CASE_PAR
#undef CHECK_INIT
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/model-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/viruses-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_VIRUSES_BONES_HPP
#define EPIWORLD_VIRUSES_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;

/**
 * @brief Set of viruses (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Viruses {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< VirusPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_viruses;

public:

    Viruses() = delete;
    Viruses(Agent<TSeq> & p) : dat(&p.viruses), n_viruses(&p.n_viruses) {};

    typename std::vector< VirusPtr<TSeq> >::iterator begin();
    typename std::vector< VirusPtr<TSeq> >::iterator end();

    VirusPtr<TSeq> & operator()(size_t i);
    VirusPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

};

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::iterator Viruses<TSeq>::begin()
{

    if (*n_viruses == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::iterator Viruses<TSeq>::end()
{
    
    #ifdef EPI_DEBUG
    if (dat->size() < *n_viruses)
        throw EPI_DEBUG_ERROR(std::logic_error, "Viruses:: The end of the virus is out of range");
    #endif 

    return begin() + *n_viruses;
}

template<typename TSeq>
inline VirusPtr<TSeq> & Viruses<TSeq>::operator()(size_t i)
{

    if (i >= *n_viruses)
        throw std::range_error("Virus index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline VirusPtr<TSeq> & Viruses<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Viruses<TSeq>::size() const noexcept 
{
    return *n_viruses;
}

/**
 * @brief Set of Viruses (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Viruses_const {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< VirusPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_viruses;

public:

    Viruses_const() = delete;
    Viruses_const(const Agent<TSeq> & p) : dat(&p.viruses), n_viruses(&p.n_viruses) {};

    typename std::vector< VirusPtr<TSeq> >::const_iterator begin() const;
    typename std::vector< VirusPtr<TSeq> >::const_iterator end() const;

    const VirusPtr<TSeq> & operator()(size_t i);
    const VirusPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

};

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::const_iterator Viruses_const<TSeq>::begin() const {

    if (*n_viruses == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::const_iterator Viruses_const<TSeq>::end() const {

    #ifdef EPI_DEBUG
    if (dat->size() < *n_viruses)
        throw EPI_DEBUG_ERROR(std::logic_error, "Viruses_const:: The end of the virus is out of range");
    #endif 
    return begin() + *n_viruses;
}

template<typename TSeq>
inline const VirusPtr<TSeq> & Viruses_const<TSeq>::operator()(size_t i)
{

    if (i >= *n_viruses)
        throw std::range_error("Virus index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline const VirusPtr<TSeq> & Viruses_const<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Viruses_const<TSeq>::size() const noexcept 
{
    return *n_viruses;
}



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/viruses-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/virus-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_VIRUS_HPP
#define EPIWORLD_VIRUS_HPP

template<typename TSeq>
class Agent;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Model;

/**
 * @brief Virus
 * 
 * @tparam TSeq 
 * @details
 * Raw transmisibility of a virus should be a function of its genetic
 * sequence. Nonetheless, transmisibility can be reduced as a result of
 * having one or more tools to fight the virus. Because of this, transmisibility
 * should be a function of the agent.
 */
template<typename TSeq>
class Virus {
    friend class Agent<TSeq>;
    friend class Model<TSeq>;
    friend class DataBase<TSeq>;
    friend void default_add_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:
    
    Agent<TSeq> * agent       = nullptr;
    int       pos_in_agent    = -99; ///< Location in the agent
    int agent_exposure_number = -99;

    std::shared_ptr<TSeq> baseline_sequence = nullptr;
    std::shared_ptr<std::string> virus_name = nullptr;
    int date = -99;
    int id   = -99;
    bool active = true;
    MutFun<TSeq>          mutation_fun                 = nullptr;
    PostRecoveryFun<TSeq> post_recovery_fun            = nullptr;
    VirusFun<TSeq>        probability_of_infecting_fun = nullptr;
    VirusFun<TSeq>        probability_of_recovery_fun  = nullptr;
    VirusFun<TSeq>        probability_of_death_fun     = nullptr;

    // Setup parameters
    std::vector< epiworld_double * > params;
    std::vector< epiworld_double > data;

    epiworld_fast_int status_init    = -99; ///< Change of state when added to agent.
    epiworld_fast_int status_post    = -99; ///< Change of state when removed from agent.
    epiworld_fast_int status_removed = -99; ///< Change of state when agent is removed

    epiworld_fast_int queue_init    = QueueValues::Everyone; ///< Change of state when added to agent.
    epiworld_fast_int queue_post    = -QueueValues::Everyone; ///< Change of state when removed from agent.
    epiworld_fast_int queue_removed = -99; ///< Change of state when agent is removed

public:
    Virus(std::string name = "unknown virus");

    void mutate(Model<TSeq> * model);
    void set_mutation(MutFun<TSeq> fun);
    
    const TSeq* get_sequence();
    void set_sequence(TSeq sequence);
    
    Agent<TSeq> * get_agent();
    void set_agent(Agent<TSeq> * p, epiworld_fast_uint idx);
    
    void set_date(int d);
    int get_date() const;

    void set_id(int idx);
    int get_id() const;

    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_prob_infecting(Model<TSeq> * model);
    epiworld_double get_prob_recovery(Model<TSeq> * model);
    epiworld_double get_prob_death(Model<TSeq> * model);
    
    void post_recovery(Model<TSeq> * model);
    void set_post_recovery(PostRecoveryFun<TSeq> fun);
    void set_post_immunity(epiworld_double prob);
    void set_post_immunity(epiworld_double * prob);

    void set_prob_infecting_fun(VirusFun<TSeq> fun);
    void set_prob_recovery_fun(VirusFun<TSeq> fun);
    void set_prob_death_fun(VirusFun<TSeq> fun);
    
    void set_prob_infecting(epiworld_double * prob);
    void set_prob_recovery(epiworld_double * prob);
    void set_prob_death(epiworld_double * prob);
    
    void set_prob_infecting(epiworld_double prob);
    void set_prob_recovery(epiworld_double prob);
    void set_prob_death(epiworld_double prob);
    ///@}


    void set_name(std::string name);
    std::string get_name() const;

    std::vector< epiworld_double > & get_data();

    /**
     * @name Get and set the state and queue
     * 
     * After applied, viruses can change the state and affect
     * the queue of agents. These function sets the default values,
     * which are retrieved when adding or removing a virus does not
     * specify a change in state or in queue.
     * 
     * @param init After the virus/tool is added to the agent.
     * @param end After the virus/tool is removed.
     * @param removed After the agent (Agent) is removed.
     */
    ///@{
    void set_state(
        epiworld_fast_int init,
        epiworld_fast_int end,
        epiworld_fast_int removed = -99
        );
        
    void set_queue(
        epiworld_fast_int init,
        epiworld_fast_int end,
        epiworld_fast_int removed = -99
        );

    void get_state(
        epiworld_fast_int * init,
        epiworld_fast_int * end,
        epiworld_fast_int * removed = nullptr
        );

    void get_queue(
        epiworld_fast_int * init,
        epiworld_fast_int * end,
        epiworld_fast_int * removed = nullptr
        );
    ///@}

    bool operator==(const Virus<TSeq> & other) const;
    bool operator!=(const Virus<TSeq> & other) const {return !operator==(other);};

    void print() const;

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/virus-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/virus-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_VIRUS_MEAT_HPP
#define EPIWORLD_VIRUS_MEAT_HPP

template<typename TSeq>
inline Virus<TSeq>::Virus(std::string name) {
    set_name(name);
}

// template<typename TSeq>
// inline Virus<TSeq>::Virus(TSeq sequence, std::string name) {
//     baseline_sequence = std::make_shared<TSeq>(sequence);
//     set_name(name);
// }

template<typename TSeq>
inline void Virus<TSeq>::mutate(
    Model<TSeq> * model
) {

    if (mutation_fun)
        if (mutation_fun(agent, *this, model))
            model->get_db().record_variant(*this);

    return;
    
}

template<typename TSeq>
inline void Virus<TSeq>::set_mutation(
    MutFun<TSeq> fun
) {
    mutation_fun = MutFun<TSeq>(fun);
}

template<typename TSeq>
inline const TSeq * Virus<TSeq>::get_sequence()
{

    return &(*baseline_sequence);

}

template<typename TSeq>
inline void Virus<TSeq>::set_sequence(TSeq sequence)
{

    baseline_sequence = std::make_shared<TSeq>(sequence);
    return;

}

template<typename TSeq>
inline Agent<TSeq> * Virus<TSeq>::get_agent()
{

    return agent;

}

template<typename TSeq>
inline void Virus<TSeq>::set_agent(Agent<TSeq> * p, epiworld_fast_uint idx)
{

    #ifdef EPI_DEBUG
    if (idx >= p->viruses.size())
    {
        printf_epiworld(
            "[epi-debug]Virus::set_agent id to set up is outside of range."
            );
    }
    #endif

    agent        = p;
    pos_in_agent = static_cast<int>(idx);

}

template<typename TSeq>
inline void Virus<TSeq>::set_id(int idx)
{

    id = idx;
    return;

}

template<typename TSeq>
inline int Virus<TSeq>::get_id() const
{
    
    return id;

}

template<typename TSeq>
inline void Virus<TSeq>::set_date(int d) 
{

    date = d;
    return;

}

template<typename TSeq>
inline int Virus<TSeq>::get_date() const
{
    
    return date;
    
}

template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_infecting(
    Model<TSeq> * model
)
{

    if (probability_of_infecting_fun)
        return probability_of_infecting_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_INFECTION;

}



template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_recovery(
    Model<TSeq> * model
)
{

    if (probability_of_recovery_fun)
        return probability_of_recovery_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_RECOVERY;

}



template<typename TSeq>
inline epiworld_double Virus<TSeq>::get_prob_death(
    Model<TSeq> * model
)
{

    if (probability_of_death_fun)
        return probability_of_death_fun(agent, *this, model);
        
    return EPI_DEFAULT_VIRUS_PROB_DEATH;

}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting_fun(VirusFun<TSeq> fun)
{
    probability_of_infecting_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery_fun(VirusFun<TSeq> fun)
{
    probability_of_recovery_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death_fun(VirusFun<TSeq> fun)
{
    probability_of_death_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting(epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_infecting_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery(epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_recovery_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death(epiworld_double * prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return *prob;
        };
    
    probability_of_death_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_infecting(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_infecting_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_recovery(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_recovery_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_prob_death(epiworld_double prob)
{
    VirusFun<TSeq> tmpfun = 
        [prob](Agent<TSeq> *, Virus<TSeq> &, Model<TSeq> *)
        {
            return prob;
        };
    
    probability_of_death_fun = tmpfun;
}

template<typename TSeq>
inline void Virus<TSeq>::set_post_recovery(PostRecoveryFun<TSeq> fun)
{
    if (post_recovery_fun)
    {
        printf_epiworld(
            "Warning: a PostRecoveryFun is alreay in place (overwriting)."
            );
    }

    post_recovery_fun = fun;
}

template<typename TSeq>
inline void Virus<TSeq>::post_recovery(
    Model<TSeq> * model
)
{

    if (post_recovery_fun)
        post_recovery_fun(agent, *this, model);    

    return;
        
}

template<typename TSeq>
inline void Virus<TSeq>::set_post_immunity(
    epiworld_double prob
)
{

    if (post_recovery_fun)
    {

        std::string msg =
            std::string(
                "You cannot set post immunity when a post_recovery "
                ) +
            std::string(
                "function is already in place. Redesign the post_recovery function."
                );

        throw std::logic_error(msg);
        
    }

    // To make sure that we keep registering the virus
    ToolPtr<TSeq> __no_reinfect = std::make_shared<Tool<TSeq>>(
        "Immunity (" + *virus_name + ")"
    );

    __no_reinfect->set_susceptibility_reduction(prob);
    __no_reinfect->set_death_reduction(0.0);
    __no_reinfect->set_transmission_reduction(0.0);
    __no_reinfect->set_recovery_enhancer(0.0);

    PostRecoveryFun<TSeq> tmpfun = 
        [__no_reinfect](
            Agent<TSeq> * p, Virus<TSeq> &, Model<TSeq> * m
            )
        {
            
            // Have we registered the tool?
            if (__no_reinfect->get_id() == -99)
                m->get_db().record_tool(*__no_reinfect);

            p->add_tool(__no_reinfect, m);

            return;

        };

    post_recovery_fun = tmpfun;

}

template<typename TSeq>
inline void Virus<TSeq>::set_post_immunity(
    epiworld_double * prob
)
{

    if (post_recovery_fun)
    {

        std::string msg =
            std::string(
                "You cannot set post immunity when a post_recovery "
                ) +
            std::string(
                "function is already in place. Redesign the post_recovery function."
                );

        throw std::logic_error(msg);

    }

    // To make sure that we keep registering the virus
    ToolPtr<TSeq> __no_reinfect = std::make_shared<Tool<TSeq>>(
        "Immunity (" + *virus_name + ")"
    );

    __no_reinfect->set_susceptibility_reduction(prob);
    __no_reinfect->set_death_reduction(0.0);
    __no_reinfect->set_transmission_reduction(0.0);
    __no_reinfect->set_recovery_enhancer(0.0);

    PostRecoveryFun<TSeq> tmpfun = 
        [__no_reinfect](Agent<TSeq> * p, Virus<TSeq> &, Model<TSeq> * m)
        {

            // Have we registered the tool?
            if (__no_reinfect->get_id() == -99)
                m->get_db().record_tool(*__no_reinfect);

            p->add_tool(__no_reinfect, m);

            return;

        };

    post_recovery_fun = tmpfun;

}

template<typename TSeq>
inline void Virus<TSeq>::set_name(std::string name)
{

    if (name == "")
        virus_name = nullptr;
    else
        virus_name = std::make_shared<std::string>(name);

}

template<typename TSeq>
inline std::string Virus<TSeq>::get_name() const
{

    if (virus_name)
        return *virus_name;
    
    return "unknown virus";

}

template<typename TSeq>
inline std::vector< epiworld_double > & Virus<TSeq>::get_data() {
    return data;
}

template<typename TSeq>
inline void Virus<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end,
    epiworld_fast_int removed
)
{
    status_init    = init;
    status_post    = end;
    status_removed = removed;
}

template<typename TSeq>
inline void Virus<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end,
    epiworld_fast_int removed
)
{

    queue_init    = init;
    queue_post     = end;
    queue_removed = removed;

}

template<typename TSeq>
inline void Virus<TSeq>::get_state(
    epiworld_fast_int * init,
    epiworld_fast_int * end,
    epiworld_fast_int * removed
)
{

    if (init != nullptr)
        *init = status_init;

    if (end != nullptr)
        *end = status_post;

    if (removed != nullptr)
        *removed = status_removed;

}

template<typename TSeq>
inline void Virus<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * end,
    epiworld_fast_int * removed
)
{

    if (init != nullptr)
        *init = queue_init;

    if (end != nullptr)
        *end = queue_post;

    if (removed != nullptr)
        *removed = queue_removed;
        
}

template<>
inline bool Virus<std::vector<int>>::operator==(
    const Virus<std::vector<int>> & other
    ) const
{
    

    EPI_DEBUG_FAIL_AT_TRUE(
        baseline_sequence->size() != other.baseline_sequence->size(),
        "Virus:: baseline_sequence don't match"
        )

    for (size_t i = 0u; i < baseline_sequence->size(); ++i)
    {

        EPI_DEBUG_FAIL_AT_TRUE(
            baseline_sequence->operator[](i) != other.baseline_sequence->operator[](i),
            "Virus:: baseline_sequence[i] don't match"
            )

    }

    EPI_DEBUG_FAIL_AT_TRUE(
        virus_name != other.virus_name,
        "Virus:: virus_name don't match"
        )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        status_init != other.status_init,
        "Virus:: status_init don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        status_post != other.status_post,
        "Virus:: status_post don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        status_removed != other.status_removed,
        "Virus:: status_removed don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_init != other.queue_init,
        "Virus:: queue_init don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_post != other.queue_post,
        "Virus:: queue_post don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_removed != other.queue_removed,
        "Virus:: queue_removed don't match"
        )

    return true;

}

template<typename TSeq>
inline bool Virus<TSeq>::operator==(const Virus<TSeq> & other) const
{
    
    EPI_DEBUG_FAIL_AT_TRUE(
        *baseline_sequence != *other.baseline_sequence,
        "Virus:: baseline_sequence don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        virus_name != other.virus_name,
        "Virus:: virus_name don't match"
    )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        status_init != other.status_init,
        "Virus:: status_init don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        status_post != other.status_post,
        "Virus:: status_post don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        status_removed != other.status_removed,
        "Virus:: status_removed don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_init != other.queue_init,
        "Virus:: queue_init don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_post != other.queue_post,
        "Virus:: queue_post don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        queue_removed != other.queue_removed,
        "Virus:: queue_removed don't match"
    )

    return true;

}

template<typename TSeq>
inline void Virus<TSeq>::print() const
{

    printf_epiworld("Virus          : %s\n", virus_name->c_str());
    printf_epiworld("status_init    : %i\n", status_init);
    printf_epiworld("status_post    : %i\n", status_post);
    printf_epiworld("status_removed : %i\n", status_removed);
    printf_epiworld("queue_init     : %i\n", queue_init);
    printf_epiworld("queue_post     : %i\n", queue_post);
    printf_epiworld("queue_removed  : %i\n", queue_removed);

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/virus-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


    
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/tools-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_TOOLS_BONES_HPP
#define EPIWORLD_TOOLS_BONES_HPP

template<typename TSeq>
class Tool;

template<typename TSeq>
class Agent;

// #define ToolPtr<TSeq> std::shared_ptr< Tool<TSeq> >

/**
 * @brief Set of tools (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Tools {
    friend class Tool<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< ToolPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_tools;

public:

    Tools() = delete;
    Tools(Agent<TSeq> & p) : dat(&p.tools), n_tools(&p.n_tools) {};

    typename std::vector< ToolPtr<TSeq> >::iterator begin();
    typename std::vector< ToolPtr<TSeq> >::iterator end();

    ToolPtr<TSeq> & operator()(size_t i);
    ToolPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

};

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::iterator Tools<TSeq>::begin()
{

    if (*n_tools == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::iterator Tools<TSeq>::end()
{
     
    return begin() + *n_tools;
}

template<typename TSeq>
inline ToolPtr<TSeq> & Tools<TSeq>::operator()(size_t i)
{

    if (i >= *n_tools)
        throw std::range_error("Tool index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline ToolPtr<TSeq> & Tools<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Tools<TSeq>::size() const noexcept 
{
    return *n_tools;
}

/**
 * @brief Set of Tools (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Tools_const {
    friend class Tool<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< ToolPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_tools;

public:

    Tools_const() = delete;
    Tools_const(const Agent<TSeq> & p) : dat(&p.tools), n_tools(&p.n_tools) {};

    typename std::vector< ToolPtr<TSeq> >::const_iterator begin() const;
    typename std::vector< ToolPtr<TSeq> >::const_iterator end() const;

    const ToolPtr<TSeq> & operator()(size_t i);
    const ToolPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

};

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::const_iterator Tools_const<TSeq>::begin() const {

    if (*n_tools == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::const_iterator Tools_const<TSeq>::end() const {
     
    return begin() + *n_tools;
}

template<typename TSeq>
inline const ToolPtr<TSeq> & Tools_const<TSeq>::operator()(size_t i)
{

    if (i >= *n_tools)
        throw std::range_error("Tool index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline const ToolPtr<TSeq> & Tools_const<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Tools_const<TSeq>::size() const noexcept 
{
    return *n_tools;
}



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/tools-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/tool-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



#ifndef EPIWORLD_TOOL_BONES_HPP
#define EPIWORLD_TOOL_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;

template<typename TSeq>
class Model;

template<typename TSeq>
class Tool;

/**
 * @brief Tools for defending the agent against the virus
 * 
 * @tparam TSeq Type of sequence
 */
template<typename TSeq> 
class Tool {
    friend class Agent<TSeq>;
    friend class Model<TSeq>;
    friend void default_add_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:

    Agent<TSeq> * agent = nullptr;
    int pos_in_agent        = -99; ///< Location in the agent

    int date = -99;
    int id   = -99;
    std::shared_ptr<std::string> tool_name     = nullptr;
    std::shared_ptr<TSeq> sequence             = nullptr;
    ToolFun<TSeq> susceptibility_reduction_fun = nullptr;
    ToolFun<TSeq> transmission_reduction_fun   = nullptr;
    ToolFun<TSeq> recovery_enhancer_fun        = nullptr;
    ToolFun<TSeq> death_reduction_fun          = nullptr;

    // Setup parameters
    std::vector< epiworld_double * > params;  

    epiworld_fast_int status_init = -99;
    epiworld_fast_int status_post = -99;

    epiworld_fast_int queue_init = QueueValues::NoOne; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = QueueValues::NoOne; ///< Change of state when removed from agent.

    void set_agent(Agent<TSeq> * p, size_t idx);

public:
    Tool(std::string name = "unknown tool");
    // Tool(TSeq d, std::string name = "unknown tool");

    void set_sequence(TSeq d);
    void set_sequence(std::shared_ptr<TSeq> d);
    std::shared_ptr<TSeq> get_sequence();

    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_susceptibility_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_transmission_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_recovery_enhancer(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_death_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    
    void set_susceptibility_reduction_fun(ToolFun<TSeq> fun);
    void set_transmission_reduction_fun(ToolFun<TSeq> fun);
    void set_recovery_enhancer_fun(ToolFun<TSeq> fun);
    void set_death_reduction_fun(ToolFun<TSeq> fun);

    void set_susceptibility_reduction(epiworld_double * prob);
    void set_transmission_reduction(epiworld_double * prob);
    void set_recovery_enhancer(epiworld_double * prob);
    void set_death_reduction(epiworld_double * prob);

    void set_susceptibility_reduction(epiworld_double prob);
    void set_transmission_reduction(epiworld_double prob);
    void set_recovery_enhancer(epiworld_double prob);
    void set_death_reduction(epiworld_double prob);
    ///@}

    void set_name(std::string name);
    std::string get_name() const;

    Agent<TSeq> * get_agent();
    int get_id() const;
    void set_id(int id);
    void set_date(int d);
    int get_date() const;

    void set_state(epiworld_fast_int init, epiworld_fast_int post);
    void set_queue(epiworld_fast_int init, epiworld_fast_int post);
    void get_state(epiworld_fast_int * init, epiworld_fast_int * post);
    void get_queue(epiworld_fast_int * init, epiworld_fast_int * post);

    bool operator==(const Tool<TSeq> & other) const;
    bool operator!=(const Tool<TSeq> & other) const {return !operator==(other);};

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/tool-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/tool-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



#ifndef EPIWORLD_TOOLS_MEAT_HPP
#define EPIWORLD_TOOLS_MEAT_HPP

template<typename TSeq>
inline Tool<TSeq>::Tool(std::string name)
{
    set_name(name);
}

// template<typename TSeq>
// inline Tool<TSeq>::Tool(TSeq d, std::string name) {
//     sequence = std::make_shared<TSeq>(d);
//     tool_name = std::make_shared<std::string>(name);
// }

template<typename TSeq>
inline void Tool<TSeq>::set_sequence(TSeq d) {
    sequence = std::make_shared<TSeq>(d);
}

template<typename TSeq>
inline void Tool<TSeq>::set_sequence(std::shared_ptr<TSeq> d) {
    sequence = d;
}

template<typename TSeq>
inline std::shared_ptr<TSeq> Tool<TSeq>::get_sequence() {
    return sequence;
}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (susceptibility_reduction_fun)
        return susceptibility_reduction_fun(*this, this->agent, v, model);

    return DEFAULT_TOOL_CONTAGION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (transmission_reduction_fun)
        return transmission_reduction_fun(*this, this->agent, v, model);

    return DEFAULT_TOOL_TRANSMISSION_REDUCTION;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (recovery_enhancer_fun)
        return recovery_enhancer_fun(*this, this->agent, v, model);

    return DEFAULT_TOOL_RECOVERY_ENHANCER;

}

template<typename TSeq>
inline epiworld_double Tool<TSeq>::get_death_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
)
{

    if (death_reduction_fun)
        return death_reduction_fun(*this, this->agent, v, model);

    return DEFAULT_TOOL_DEATH_REDUCTION;

}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction_fun(
    ToolFun<TSeq> fun
)
{
    susceptibility_reduction_fun = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction_fun(
    ToolFun<TSeq> fun
)
{
    transmission_reduction_fun = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer_fun(
    ToolFun<TSeq> fun
)
{
    recovery_enhancer_fun = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction_fun(
    ToolFun<TSeq> fun
)
{
    death_reduction_fun = fun;
}

template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    susceptibility_reduction_fun = tmpfun;

}

// EPIWORLD_SET_LAMBDA(susceptibility_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(epiworld_double * prob)
{
    
    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    transmission_reduction_fun = tmpfun;

}

// EPIWORLD_SET_LAMBDA(transmission_reduction)
template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    recovery_enhancer_fun = tmpfun;

}

// EPIWORLD_SET_LAMBDA(recovery_enhancer)
template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(epiworld_double * prob)
{

    ToolFun<TSeq> tmpfun =
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return *prob;
        };

    death_reduction_fun = tmpfun;

}

// EPIWORLD_SET_LAMBDA(death_reduction)

// #undef EPIWORLD_SET_LAMBDA
template<typename TSeq>
inline void Tool<TSeq>::set_susceptibility_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    susceptibility_reduction_fun = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_transmission_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    transmission_reduction_fun = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_recovery_enhancer(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    recovery_enhancer_fun = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_death_reduction(
    epiworld_double prob
)
{

    ToolFun<TSeq> tmpfun = 
        [prob](Tool<TSeq> &, Agent<TSeq> *, VirusPtr<TSeq>, Model<TSeq> *)
        {
            return prob;
        };

    death_reduction_fun = tmpfun;

}

template<typename TSeq>
inline void Tool<TSeq>::set_name(std::string name)
{
    if (name != "")
        tool_name = std::make_shared<std::string>(name);
}

template<typename TSeq>
inline std::string Tool<TSeq>::get_name() const {

    if (tool_name)
        return *tool_name;

    return "unknown tool";

}

template<typename TSeq>
inline Agent<TSeq> * Tool<TSeq>::get_agent()
{
    return this->agent;
}

template<typename TSeq>
inline void Tool<TSeq>::set_agent(Agent<TSeq> * p, size_t idx)
{
    agent        = p;
    pos_in_agent = static_cast<int>(idx);
}

template<typename TSeq>
inline int Tool<TSeq>::get_id() const {
    return id;
}


template<typename TSeq>
inline void Tool<TSeq>::set_id(int id)
{
    this->id = id;
}

template<typename TSeq>
inline void Tool<TSeq>::set_date(int d)
{
    this->date = d;
}

template<typename TSeq>
inline int Tool<TSeq>::get_date() const
{
    return date;
}

template<typename TSeq>
inline void Tool<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    status_init = init;
    status_post = end;
}

template<typename TSeq>
inline void Tool<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    queue_init = init;
    queue_post = end;
}

template<typename TSeq>
inline void Tool<TSeq>::get_state(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = status_init;

    if (post != nullptr)
        *post = status_post;

}

template<typename TSeq>
inline void Tool<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = queue_init;

    if (post != nullptr)
        *post = queue_post;

}

template<>
inline bool Tool<std::vector<int>>::operator==(
    const Tool<std::vector<int>> & other
    ) const
{
    
    if (sequence->size() != other.sequence->size())
        return false;

    for (size_t i = 0u; i < sequence->size(); ++i)
    {
        if (sequence->operator[](i) != other.sequence->operator[](i))
            return false;
    }

    if (tool_name != other.tool_name)
        return false;
    
    if (status_init != other.status_init)
        return false;

    if (status_post != other.status_post)
        return false;

    if (queue_init != other.queue_init)
        return false;

    if (queue_post != other.queue_post)
        return false;


    return true;

}

template<typename TSeq>
inline bool Tool<TSeq>::operator==(const Tool<TSeq> & other) const
{
    if (*sequence != *other.sequence)
        return false;

    if (tool_name != other.tool_name)
        return false;
    
    if (status_init != other.status_init)
        return false;

    if (status_post != other.status_post)
        return false;

    if (queue_init != other.queue_init)
        return false;

    if (queue_post != other.queue_post)
        return false;

    return true;

}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/tool-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/entity-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_ENTITY_BONES_HPP
#define EPIWORLD_ENTITY_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Agent;

template<typename TSeq>
class AgentsSample;

template<typename TSeq>
inline void default_add_entity(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_entity(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
class Entity {
    friend class Agent<TSeq>;
    friend class AgentsSample<TSeq>;
    friend class Model<TSeq>;
    friend void default_add_entity<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_entity<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:

    Model<TSeq> * model;
    
    int id = -1;
    std::vector< size_t > agents;   ///< Vector of agents
    std::vector< size_t > agents_location; ///< Location where the entity is stored in the agent
    size_t n_agents = 0u;

    /**
     * @name Auxiliary variables for AgentsSample<TSeq> iterators
     * 
     * @details These variables+objects are used by the AgentsSample<TSeq>
     * class for building efficient iterators over agents. The idea is to
     * reduce the memory allocation, so only during the first call of
     * AgentsSample<TSeq>::AgentsSample(Entity<TSeq>) these vectors are allocated.
     */
    ///@{
    std::vector< Agent<TSeq> * > sampled_agents;
    size_t sampled_agents_n = 0u;
    std::vector< size_t > sampled_agents_left;
    size_t sampled_agents_left_n = 0u;
    // int date_last_add_or_remove = -99; ///< Last time the entity added or removed an agent
    ///@}

    int max_capacity = -1;
    std::string entity_name = "Unknown entity";

    std::vector< epiworld_double > location = {0.0}; ///< An arbitrary vector for location

    epiworld_fast_int status_init = -99;
    epiworld_fast_int status_post = -99;

    epiworld_fast_int queue_init = 0; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = 0; ///< Change of state when removed from agent.

public:

    // Entity() = delete;
    // Entity(Entity<TSeq> & e) = delete;
    // Entity(const Entity<TSeq> & e);
    // Entity(Entity && e);
    Entity(std::string name) : entity_name(name) {};
    // Entity<TSeq> & operator=(const Entity<TSeq> & e);

    void add_agent(Agent<TSeq> & p, Model<TSeq> * model);
    void add_agent(Agent<TSeq> * p, Model<TSeq> * model);
    void rm_agent(size_t idx);
    size_t size() const noexcept;
    void set_location(std::vector< epiworld_double > loc);
    std::vector< epiworld_double > & get_location();

    typename std::vector< Agent<TSeq> * >::iterator begin();
    typename std::vector< Agent<TSeq> * >::iterator end();

    typename std::vector< Agent<TSeq> * >::const_iterator begin() const;
    typename std::vector< Agent<TSeq> * >::const_iterator end() const;

    Agent<TSeq> * operator[](size_t i);

    int get_id() const noexcept;
    const std::string & get_name() const noexcept;

    void set_state(epiworld_fast_int init, epiworld_fast_int post);
    void set_queue(epiworld_fast_int init, epiworld_fast_int post);
    void get_state(epiworld_fast_int * init, epiworld_fast_int * post);
    void get_queue(epiworld_fast_int * init, epiworld_fast_int * post);

    void reset();

    bool operator==(const Entity<TSeq> & other) const;
    bool operator!=(const Entity<TSeq> & other) const {return !operator==(other);};

};


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/entity-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/entity-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_ENTITY_MEAT_HPP
#define EPIWORLD_ENTITY_MEAT_HPP

// template<typename TSeq>
// inline Entity<TSeq>::Entity(const Entity<TSeq> & e) :
//     model(e.model),
//     id(e.id),
//     agents(0u),
//     agents_location(0u),
//     n_agents(0),
//     sampled_agents(0u),
//     sampled_agents_n(0u),
//     sampled_agents_left(0u),
//     sampled_agents_left_n(0u),
//     max_capacity(e.max_capacity),
//     entity_name(e.entity_name),
//     location(e.location),
//     status_init(e.status_init),
//     status_post(e.status_post),
//     queue_init(e.queue_init),
//     queue_post(e.queue_post)
// {

// }

template<typename TSeq>
inline void Entity<TSeq>::add_agent(
    Agent<TSeq> & p,
    Model<TSeq> * model
    )
{

    // Need to add it to the actions, through the individual
    p.add_entity(*this, model);    

}

template<typename TSeq>
inline void Entity<TSeq>::add_agent(
    Agent<TSeq> * p,
    Model<TSeq> * model
    )
{
    p->add_entity(*this, model);
}

template<typename TSeq>
inline void Entity<TSeq>::rm_agent(size_t idx)
{
    if (idx >= n_agents)
        throw std::out_of_range(
            "Trying to remove agent "+ std::to_string(idx) +
            " out of " + std::to_string(n_agents)
            );

    model->population[idx].rm_entity(*this);

    return;
}

template<typename TSeq>
inline size_t Entity<TSeq>::size() const noexcept
{
    return n_agents;
}

template<typename TSeq>
inline void Entity<TSeq>::set_location(std::vector< epiworld_double > loc)
{
    location = loc;
}

template<typename TSeq>
inline std::vector< epiworld_double > & Entity<TSeq>::get_location()
{
    return location;
}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator Entity<TSeq>::begin()
{

    if (n_agents == 0)
        return agents.end();

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator Entity<TSeq>::end()
{
    return agents.begin() + n_agents;
}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::const_iterator Entity<TSeq>::begin() const
{

    if (n_agents == 0)
        return agents.end();

    return agents.begin();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::const_iterator Entity<TSeq>::end() const
{
    return agents.begin() + n_agents;
}

template<typename TSeq>
inline Agent<TSeq> * Entity<TSeq>::operator[](size_t i)
{
    if (n_agents <= i)
        throw std::logic_error("There are not that many agents in this entity.");

    return &model->get_agents()[i];
}

template<typename TSeq>
inline int Entity<TSeq>::get_id() const noexcept
{
    return id;
}

template<typename TSeq>
inline const std::string & Entity<TSeq>::get_name() const noexcept
{
    return entity_name;
}

template<typename TSeq>
inline void Entity<TSeq>::set_state(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    status_init = init;
    status_post = end;
}

template<typename TSeq>
inline void Entity<TSeq>::set_queue(
    epiworld_fast_int init,
    epiworld_fast_int end
)
{
    queue_init = init;
    queue_post = end;
}

template<typename TSeq>
inline void Entity<TSeq>::get_state(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = status_init;

    if (post != nullptr)
        *post = status_post;

}

template<typename TSeq>
inline void Entity<TSeq>::get_queue(
    epiworld_fast_int * init,
    epiworld_fast_int * post
)
{
    if (init != nullptr)
        *init = queue_init;

    if (post != nullptr)
        *post = queue_post;

}

template<typename TSeq>
inline void Entity<TSeq>::reset()
{
    sampled_agents.clear();
    sampled_agents_n = 0u;
    sampled_agents_left.clear();
    sampled_agents_left_n = 0u;
}

template<typename TSeq>
inline bool Entity<TSeq>::operator==(const Entity<TSeq> & other) const
{

    if (id != other.id)
        return false;

    if (n_agents != other.n_agents)
        return false;

    for (size_t i = 0u; i < n_agents; ++i)
    {
        if (agents[i] != other.agents[i])
            return false;
    }


    if (max_capacity != other.max_capacity)
        return false;

    if (entity_name != other.entity_name)
        return false;

    if (location.size() != other.location.size())
        return false;

    for (size_t i = 0u; i < location.size(); ++i)
    {

        if (location[i] != other.location[i])
            return false;

    }

    if (status_init != other.status_init)
        return false;

    if (status_post != other.status_post)
        return false;

    if (queue_init != other.queue_init)
        return false;

    if (queue_post != other.queue_post)
        return false;

    return true;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/entity-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/entities-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_ENTITIES_BONES_HPP
#define EPIWORLD_ENTITIES_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;


/**
 * @brief Set of Entities (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Entities {
    friend class Entity<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< Entity<TSeq> * >  dat;
    const size_t n_entities;

public:

    Entities() = delete;
    Entities(Agent<TSeq> & p);

    typename std::vector< Entity<TSeq> * >::iterator begin();
    typename std::vector< Entity<TSeq> * >::iterator end();

    Entity<TSeq> & operator()(size_t i);
    Entity<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    bool operator==(const Entities<TSeq> & other) const;

};

template<typename TSeq>
inline Entities<TSeq>::Entities(Agent<TSeq> & p) :
    n_entities(p.get_n_entities())
{

    dat.reserve(n_entities);
    for (size_t i = 0u; i < n_entities; ++i)
        dat.push_back(&p.get_entity(i));

}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::iterator Entities<TSeq>::begin()
{

    if (n_entities == 0u)
        return dat.end();
    
    return dat.begin();
}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::iterator Entities<TSeq>::end()
{
     
    return begin() + n_entities;
}

template<typename TSeq>
inline Entity<TSeq> & Entities<TSeq>::operator()(size_t i)
{

    if (i >= n_entities)
        throw std::range_error("Entity index out of range.");

    return *dat[i];

}

template<typename TSeq>
inline Entity<TSeq> & Entities<TSeq>::operator[](size_t i)
{

    return *dat[i];

}

template<typename TSeq>
inline size_t Entities<TSeq>::size() const noexcept 
{
    return n_entities;
}

template<typename TSeq>
inline bool Entities<TSeq>::operator==(const Entities<TSeq> & other) const
{

    if (n_entities != other.n_entities)
        return false;

    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (dat[i] != other.dat[i])
            return false;
    }

    return true;
}

/**
 * @brief Set of Entities (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Entities_const {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< Entity<TSeq>* > dat;
    const size_t n_entities;

public:

    Entities_const() = delete;
    Entities_const(const Agent<TSeq> & p);

    typename std::vector< Entity<TSeq>* >::const_iterator begin();
    typename std::vector< Entity<TSeq>* >::const_iterator end();

    const Entity<TSeq> & operator()(size_t i);
    const Entity<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    bool operator==(const Entities_const<TSeq> & other) const;

};

template<typename TSeq>
inline Entities_const<TSeq>::Entities_const(const Agent<TSeq> & p) :
    n_entities(p.get_n_entities())
{

    dat.reserve(n_entities);
    for (size_t i = 0u; i < n_entities; ++i)
        dat.push_back(&p.get_entity(i));

}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::const_iterator Entities_const<TSeq>::begin() {

    if (n_entities == 0u)
        return dat.end();
    
    return dat.begin();
}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::const_iterator Entities_const<TSeq>::end() {
     
    return begin() + n_entities;
}

template<typename TSeq>
inline const Entity<TSeq> & Entities_const<TSeq>::operator()(size_t i)
{

    if (i >= n_entities)
        throw std::range_error("Entity index out of range.");

    return *dat[i];

}

template<typename TSeq>
inline const Entity<TSeq> & Entities_const<TSeq>::operator[](size_t i)
{

    return *dat[i];

}

template<typename TSeq>
inline size_t Entities_const<TSeq>::size() const noexcept 
{
    return n_entities;
}

template<typename TSeq>
inline bool Entities_const<TSeq>::operator==(const Entities_const<TSeq> & other) const
{
    
    if (n_entities != other.n_entities)
        return false;

    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (dat[i] != other.dat[i])
            return false;
    }

    return true;
}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/entities-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


    
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/agent-meat-state.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_PERSON_MEAT_STATE_HPP
#define EPIWORLD_PERSON_MEAT_STATE_HPP

// template<typename TSeq>
// class Model;

// template<typename TSeq>
// class Agent;


/**
 * @file agent-meat-state.hpp
 * @author George G. Vega Yon (g.vegayon en gmail)
 * @brief Sampling functions are getting big, so we keep them in a separate file.
 * @version 0.1
 * @date 2022-06-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//agent-meat-virus-sampling.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING
#define EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING

/**
 * @brief Functions for sampling viruses
 * 
 */
namespace sampler {

/**
 * @brief Make a function to sample from neighbors
 * 
 * This is akin to the function default_update_susceptible, with the difference
 * that it will create a function that supports excluding states from the sampling
 * frame. For example, individuals who have acquired a virus can be excluded if
 * in incubation state.
 * 
 * @tparam TSeq 
 * @param exclude unsigned vector of states that need to be excluded from the sampling
 * @return Virus<TSeq>* of the selected virus. If none selected (or none
 * available,) returns a nullptr;
 */
template<typename TSeq>
inline std::function<void(Agent<TSeq>*,Model<TSeq>*)> make_update_susceptible(
    std::vector< epiworld_fast_uint > exclude = {}
    )
{
  

    if (exclude.size() == 0u)
    {

        std::function<void(Agent<TSeq>*,Model<TSeq>*)> sampler =
            [](Agent<TSeq> * p, Model<TSeq> * m) -> void
            {

                if (p->get_n_viruses() > 0u)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has ") + std::to_string(p->get_n_viruses()) +
                        std::string(" viruses.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nvariants_tmp = 0u;
                for (auto & neighbor: p->get_neighbors()) 
                {
                            
                    for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
                    { 

                        #ifdef EPI_DEBUG
                        if (nvariants_tmp >= m->array_virus_tmp.size())
                            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                        #endif
                            
                        /* And it is a function of susceptibility_reduction as well */ 
                        m->array_double_tmp[nvariants_tmp] =
                            (1.0 - p->get_susceptibility_reduction(v, m)) * 
                            v->get_prob_infecting(m) * 
                            (1.0 - neighbor->get_transmission_reduction(v, m)) 
                            ; 
                    
                        m->array_virus_tmp[nvariants_tmp++] = &(*v);
                        
                    } 
                }

                // No virus to compute
                if (nvariants_tmp == 0u)
                    return;

                // Running the roulette
                int which = roulette(nvariants_tmp, m);

                if (which < 0)
                    return;

                p->add_virus(*m->array_virus_tmp[which], m);

                return; 
            };

        return sampler;

    } else {

        // Making room for the query
        std::shared_ptr<std::vector<bool>> exclude_agent_bool =
            std::make_shared<std::vector<bool>>(0);

        std::shared_ptr<std::vector<epiworld_fast_uint>> exclude_agent_bool_idx =
            std::make_shared<std::vector<epiworld_fast_uint>>(exclude);

        std::function<void(Agent<TSeq>*,Model<TSeq>*)> sampler =
            [exclude_agent_bool,exclude_agent_bool_idx](Agent<TSeq> * p, Model<TSeq> * m) -> void
            {

                // The first time we call it, we need to initialize the vector
                if (exclude_agent_bool->size() == 0u)
                {

                    exclude_agent_bool->resize(m->get_state().size(), false);
                    for (auto s : *exclude_agent_bool_idx)
                    {
                        if (s >= exclude_agent_bool->size())
                            throw std::logic_error(
                                std::string("You are trying to exclude a state that is out of range: ") +
                                std::to_string(s) + std::string(". There are only ") +
                                std::to_string(exclude_agent_bool->size()) + 
                                std::string(" states in the model.")
                                );

                        exclude_agent_bool->operator[](s) = true;

                    }

                }                    

                if (p->get_n_viruses() > 0u)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has ") + std::to_string(p->get_n_viruses()) +
                        std::string(" viruses.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nvariants_tmp = 0u;
                for (auto & neighbor: p->get_neighbors()) 
                {

                    // If the state is in the list, exclude it
                    if (exclude_agent_bool->operator[](neighbor->get_state()))
                        continue;
                            
                    for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
                    { 

                        #ifdef EPI_DEBUG
                        if (nvariants_tmp >= m->array_virus_tmp.size())
                            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                            
                        #endif
                            
                        /* And it is a function of susceptibility_reduction as well */ 
                        m->array_double_tmp[nvariants_tmp] =
                            (1.0 - p->get_susceptibility_reduction(v, m)) * 
                            v->get_prob_infecting(m) * 
                            (1.0 - neighbor->get_transmission_reduction(v, m)) 
                            ; 
                    
                        m->array_virus_tmp[nvariants_tmp++] = &(*v);
                        
                    } 
                }

                // No virus to compute
                if (nvariants_tmp == 0u)
                    return;

                // Running the roulette
                int which = roulette(nvariants_tmp, m);

                if (which < 0)
                    return;

                p->add_virus(*m->array_virus_tmp[which], m); 

                return;

            };

        return sampler;

    }
    
}

/**
 * @brief Make a function to sample from neighbors
 * 
 * This is akin to the function default_update_susceptible, with the difference
 * that it will create a function that supports excluding states from the sampling
 * frame. For example, individuals who have acquired a virus can be excluded if
 * in incubation state.
 * 
 * @tparam TSeq 
 * @param exclude unsigned vector of states that need to be excluded from the sampling
 * @return Virus<TSeq>* of the selected virus. If none selected (or none
 * available,) returns a nullptr;
 */
template<typename TSeq = int>
inline std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> make_sample_virus_neighbors(
    std::vector< epiworld_fast_uint > exclude = {}
)
{
    if (exclude.size() == 0u)
    {

        std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> res = 
            [](Agent<TSeq> * p, Model<TSeq> * m) -> Virus<TSeq>* {

                if (p->get_n_viruses() > 0u)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has ") + std::to_string(p->get_n_viruses()) +
                        std::string(" viruses.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nvariants_tmp = 0u;
                for (auto & neighbor: p->get_neighbors()) 
                {
                            
                    for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
                    { 

                        #ifdef EPI_DEBUG
                        if (nvariants_tmp >= m->array_virus_tmp.size())
                            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                        #endif
                            
                        /* And it is a function of susceptibility_reduction as well */ 
                        m->array_double_tmp[nvariants_tmp] =
                            (1.0 - p->get_susceptibility_reduction(v, m)) * 
                            v->get_prob_infecting(m) * 
                            (1.0 - neighbor->get_transmission_reduction(v, m)) 
                            ; 
                    
                        m->array_virus_tmp[nvariants_tmp++] = &(*v);
                        
                    } 
                }

                // No virus to compute
                if (nvariants_tmp == 0u)
                    return nullptr;

                // Running the roulette
                int which = roulette(nvariants_tmp, m);

                if (which < 0)
                    return nullptr;

                return m->array_virus_tmp[which]; 

            };

        return res;


    } else {

        // Making room for the query
        std::shared_ptr<std::vector<bool>> exclude_agent_bool =
            std::make_shared<std::vector<bool>>(0);

        std::shared_ptr<std::vector<epiworld_fast_uint>> exclude_agent_bool_idx =
            std::make_shared<std::vector<epiworld_fast_uint>>(exclude);


        std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> res = 
            [exclude_agent_bool,exclude_agent_bool_idx](Agent<TSeq> * p, Model<TSeq> * m) -> Virus<TSeq>* {

                // The first time we call it, we need to initialize the vector
                if (exclude_agent_bool->size() == 0u)
                {

                    exclude_agent_bool->resize(m->get_state().size(), false);
                    for (auto s : *exclude_agent_bool_idx)
                    {
                        if (s >= exclude_agent_bool->size())
                            throw std::logic_error(
                                std::string("You are trying to exclude a state that is out of range: ") +
                                std::to_string(s) + std::string(". There are only ") +
                                std::to_string(exclude_agent_bool->size()) + 
                                std::string(" states in the model.")
                                );

                        exclude_agent_bool->operator[](s) = true;

                    }

                }    
                
                if (p->get_n_viruses() > 0u)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has ") + std::to_string(p->get_n_viruses()) +
                        std::string(" viruses.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nvariants_tmp = 0u;
                for (auto & neighbor: p->get_neighbors()) 
                {

                    // If the state is in the list, exclude it
                    if (exclude_agent_bool->operator[](neighbor->get_state()))
                        continue;
                            
                    for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
                    { 

                        #ifdef EPI_DEBUG
                        if (nvariants_tmp >= m->array_virus_tmp.size())
                            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                        #endif
                            
                        /* And it is a function of susceptibility_reduction as well */ 
                        m->array_double_tmp[nvariants_tmp] =
                            (1.0 - p->get_susceptibility_reduction(v, m)) * 
                            v->get_prob_infecting(m) * 
                            (1.0 - neighbor->get_transmission_reduction(v, m)) 
                            ; 
                    
                        m->array_virus_tmp[nvariants_tmp++] = &(*v);
                        
                    } 
                }

                // No virus to compute
                if (nvariants_tmp == 0u)
                    return nullptr;

                // Running the roulette
                int which = roulette(nvariants_tmp, m);

                if (which < 0)
                    return nullptr;

                return m->array_virus_tmp[which]; 

            };

        return res;

    }

}

/**
 * @brief Sample from neighbors pool of viruses (at most one)
 * 
 * This function samples at most one virus from the pool of
 * viruses from its neighbors. If no virus is selected, the function
 * returns a `nullptr`, otherwise it returns a pointer to the
 * selected virus.
 * 
 * This can be used to build a new update function (EPI_NEW_UPDATEFUN.)
 * 
 * @tparam TSeq 
 * @param p Pointer to person 
 * @param m Pointer to the model
 * @return Virus<TSeq>* of the selected virus. If none selected (or none
 * available,) returns a nullptr;
 */
template<typename TSeq = int>
inline Virus<TSeq> * sample_virus_single(Agent<TSeq> * p, Model<TSeq> * m)
{

    if (p->get_n_viruses() > 0u)
        throw std::logic_error(
            std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
            std::string("Agent id ") + std::to_string(p->get_id()) +
            std::string(" has ") + std::to_string(p->get_n_viruses()) +
            std::string(" viruses.")
            );

    // This computes the prob of getting any neighbor variant
    size_t nvariants_tmp = 0u;
    for (auto & neighbor: p->get_neighbors()) 
    {   
        #ifdef EPI_DEBUG
        int _vcount_neigh = 0;
        #endif                 
        for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
        { 

            #ifdef EPI_DEBUG
            if (nvariants_tmp >= m->array_virus_tmp.size())
                throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
            #endif
                
            /* And it is a function of susceptibility_reduction as well */ 
            m->array_double_tmp[nvariants_tmp] =
                (1.0 - p->get_susceptibility_reduction(v, m)) * 
                v->get_prob_infecting(m) * 
                (1.0 - neighbor->get_transmission_reduction(v, m)) 
                ; 
        
            m->array_virus_tmp[nvariants_tmp++] = &(*v);

            #ifdef EPI_DEBUG
            if (
                (m->array_double_tmp[nvariants_tmp - 1] < 0.0) |
                (m->array_double_tmp[nvariants_tmp - 1] > 1.0)
                )
            {
                printf_epiworld(
                    "[epi-debug] Agent %i's virus %i has transmission prob outside of [0, 1]: %.4f!\n",
                    static_cast<int>(neighbor->get_id()),
                    static_cast<int>(_vcount_neigh++),
                    m->array_double_tmp[nvariants_tmp - 1]
                    );
            }
            #endif
            
        } 
    }


    // No virus to compute
    if (nvariants_tmp == 0u)
        return nullptr;

    #ifdef EPI_DEBUG
    m->get_db().n_transmissions_potential++;
    #endif

    // Running the roulette
    int which = roulette(nvariants_tmp, m);

    if (which < 0)
        return nullptr;

    #ifdef EPI_DEBUG
    m->get_db().n_transmissions_today++;
    #endif

    return m->array_virus_tmp[which]; 
    
}

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//agent-meat-virus-sampling.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/




template<typename TSeq = EPI_DEFAULT_TSEQ>
inline void default_update_susceptible(
    Agent<TSeq> * p,
    Model<TSeq> * m
    )
{

    Virus<TSeq> * virus = sampler::sample_virus_single<TSeq>(p, m);
    
    if (virus == nullptr)
        return;

    p->add_virus(*virus, m); 

    return;

}

template<typename TSeq = EPI_DEFAULT_TSEQ>
inline void default_update_exposed(Agent<TSeq> * p, Model<TSeq> * m) {

    if (p->get_n_viruses() == 0u)
        throw std::logic_error(
            std::string("Using the -default_update_exposed- on agents WITHOUT viruses makes no sense! ") +
            std::string("Agent id ") + std::to_string(p->get_id()) + std::string(" has no virus registered.")
            );

    // Odd: Die, Even: Recover
    epiworld_fast_uint n_events = 0u;
    for (const auto & v : p->get_viruses())
    {

        // Die
        m->array_double_tmp[n_events++] = 
            v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, m)); 

        // Recover
        m->array_double_tmp[n_events++] = 
            1.0 - (1.0 - v->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(v, m)); 

    }

    #ifdef EPI_DEBUG
    if (n_events == 0u)
    {
        printf_epiworld(
            "[epi-debug] agent %i has 0 possible events!!\n",
            static_cast<int>(p->get_id())
            );
        throw std::logic_error("Zero events in exposed.");
    }
    #else
    if (n_events == 0u)
        return;
    #endif
    

    // Running the roulette
    int which = roulette(n_events, m);

    if (which < 0)
        return;

    // Which roulette happen?
    if ((which % 2) == 0) // If odd
    {

        size_t which_v = std::ceil(which / 2);
        p->rm_agent_by_virus(which_v, m);
        
    } else {

        size_t which_v = std::floor(which / 2);
        p->rm_virus(which_v, m);

    }

    return ;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/agent-meat-state.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/agent-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_PERSON_BONES_HPP
#define EPIWORLD_PERSON_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Viruses;

template<typename TSeq>
class Viruses_const;

template<typename TSeq>
class Tool;

template<typename TSeq>
class Tools;

template<typename TSeq>
class Tools_const;

template<typename TSeq>
class Queue;

template<typename TSeq>
struct Action;

template<typename TSeq>
class Entity;

template<typename TSeq>
class Entities;

template<typename TSeq>
inline void default_add_virus(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_tool(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_entity(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_virus(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_tool(Action<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_entity(Action<TSeq> & a, Model<TSeq> * m);



/**
 * @brief Agent (agents)
 * 
 * @tparam TSeq Sequence type (should match `TSeq` across the model)
 */
template<typename TSeq>
class Agent {
    friend class Model<TSeq>;
    friend class Virus<TSeq>;
    friend class Viruses<TSeq>;
    friend class Viruses_const<TSeq>;
    friend class Tool<TSeq>;
    friend class Tools<TSeq>;
    friend class Tools_const<TSeq>;
    friend class Queue<TSeq>;
    friend class Entities<TSeq>;
    friend class AgentsSample<TSeq>;
    friend void default_add_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_add_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_add_entity<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_entity<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:
    
    Model<TSeq> * model;

    std::vector< size_t > neighbors;
    std::vector< size_t > neighbors_locations;
    size_t n_neighbors = 0u;

    std::vector< size_t > entities;
    std::vector< size_t > entities_locations;
    size_t n_entities = 0u;

    epiworld_fast_uint state = 0u;
    epiworld_fast_uint status_prev = 0u; ///< For accounting, if need to undo a change.
    
    int status_last_changed = -1; ///< Last time the agent was updated.
    int id = -1;
    
    std::vector< VirusPtr<TSeq> > viruses;
    epiworld_fast_uint n_viruses = 0u;

    std::vector< ToolPtr<TSeq> > tools;
    epiworld_fast_uint n_tools = 0u;

    ActionFun<TSeq> add_virus_  = default_add_virus<TSeq>;
    ActionFun<TSeq> add_tool_   = default_add_tool<TSeq>;
    ActionFun<TSeq> add_entity_ = default_add_entity<TSeq>;

    ActionFun<TSeq> rm_virus_  = default_rm_virus<TSeq>;
    ActionFun<TSeq> rm_tool_   = default_rm_tool<TSeq>;
    ActionFun<TSeq> rm_entity_ = default_rm_entity<TSeq>;
    
    epiworld_fast_uint action_counter = 0u;

    std::vector< Agent<TSeq> * > sampled_agents;
    size_t sampled_agents_n      = 0u;
    std::vector< size_t > sampled_agents_left;
    size_t sampled_agents_left_n = 0u;
    int date_last_build_sample   = -99;

public:

    Agent();
    Agent(Agent<TSeq> && p);
    Agent(const Agent<TSeq> & p);
    Agent<TSeq> & operator=(const Agent<TSeq> & other_agent);

    /**
     * @name Add/Remove Virus/Tool
     * 
     * Any of these is ultimately reflected at the end of the iteration.
     * 
     * @param tool Tool to add
     * @param virus Virus to add
     * @param status_new state after the change
     * @param queue 
     */
    ///@{
    void add_tool(
        ToolPtr<TSeq> tool,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_tool(
        Tool<TSeq> tool,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_virus(
        VirusPtr<TSeq> virus,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_virus(
        Virus<TSeq> virus,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_entity(
        Entity<TSeq> & entity,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
        );

    void rm_tool(
        epiworld_fast_uint tool_idx,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_tool(
        ToolPtr<TSeq> & tool,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_virus(
        epiworld_fast_uint virus_idx,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_virus(
        VirusPtr<TSeq> & virus,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_entity(
        epiworld_fast_uint entity_idx,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_entity(
        Entity<TSeq> & entity,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_agent_by_virus(
        epiworld_fast_uint virus_idx,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    ); ///< Agent removed by virus

    void rm_agent_by_virus(
        VirusPtr<TSeq> & virus,
        Model<TSeq> * model,
        epiworld_fast_int status_new = -99,
        epiworld_fast_int queue = -99
    ); ///< Agent removed by virus
    ///@}
    
    /**
     * @name Get the rates (multipliers) for the agent
     * 
     * @param v A pointer to a virus.
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_susceptibility_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_transmission_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_recovery_enhancer(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_death_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    ///@}

    int get_id() const; ///< Id of the individual

    VirusPtr<TSeq> & get_virus(int i);
    Viruses<TSeq> get_viruses();
    const Viruses_const<TSeq> get_viruses() const;
    size_t get_n_viruses() const noexcept;

    ToolPtr<TSeq> & get_tool(int i);
    Tools<TSeq> get_tools();
    const Tools_const<TSeq> get_tools() const;
    size_t get_n_tools() const noexcept;

    void mutate_variant();
    void add_neighbor(
        Agent<TSeq> & p,
        bool check_source = true,
        bool check_target = true
        );

    /**
     * @brief Swaps neighbors between the current agent and agent `other`
     * 
     * @param other 
     * @param n_this 
     * @param n_other 
     */
    void swap_neighbors(
        Agent<TSeq> & other,
        size_t n_this,
        size_t n_other
    );

    std::vector< Agent<TSeq> * > get_neighbors();
    size_t get_n_neighbors() const;

    void change_state(
        Model<TSeq> * model,
        epiworld_fast_uint new_state,
        epiworld_fast_int queue = 0
        );

    const epiworld_fast_uint & get_state() const;

    void reset();

    bool has_tool(epiworld_fast_uint t) const;
    bool has_tool(std::string name) const;
    bool has_virus(epiworld_fast_uint t) const;
    bool has_virus(std::string name) const;

    void print(Model<TSeq> * model, bool compressed = false) const;

    /**
     * @brief Access the j-th column of the agent
     * 
     * If an external array has been specified, then these two
     * functions can be used to access additional agent's features 
     * not included in the model.
     * 
     * The `operator[]` method is with no boundary check, whereas
     * the `operator()` method checks boundaries. The former can result
     * in a segfault.
     * 
     * 
     * @param j 
     * @return double& 
     */
    ///@{
    double & operator()(size_t j);
    double & operator[](size_t j);
    double operator()(size_t j) const;
    double operator[](size_t j) const;
    ///@}

    Entities<TSeq> get_entities();
    const Entities_const<TSeq> get_entities() const;
    const Entity<TSeq> & get_entity(size_t i) const;
    Entity<TSeq> & get_entity(size_t i);
    size_t get_n_entities() const;

    bool operator==(const Agent<TSeq> & other) const;
    bool operator!=(const Agent<TSeq> & other) const {return !operator==(other);};

};



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/agent-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/agent-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_PERSON_MEAT_HPP
#define EPIWORLD_PERSON_MEAT_HPP

// template<typename Ta>
// inline bool IN(Ta & a, std::vector< Ta > & b);

#define CHECK_COALESCE_(proposed_, virus_tool_, alt_) \
    if (static_cast<int>(proposed_) == -99) {\
        if (static_cast<int>(virus_tool_) == -99) \
            (proposed_) = (alt_);\
        else (proposed_) = (virus_tool_);}

// To large to add directly here
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//agent-actions-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_AGENT_ACTIONS_MEAT_HPP
#define EPIWORLD_AGENT_ACTIONS_MEAT_HPP

template<typename TSeq>
inline void default_add_virus(Action<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> *  p = a.agent;
    VirusPtr<TSeq> v = a.virus;

    CHECK_COALESCE_(a.new_state, v->status_init, p->get_state())
    CHECK_COALESCE_(a.queue, v->queue_init, 1)

    // Has a agent? If so, we need to register the transmission
    if (v->get_agent())
    {

        // ... only if not the same agent
        if (v->get_agent()->get_id() != p->get_id())
            m->get_db().record_transmission(
                v->get_agent()->get_id(),
                p->get_id(),
                v->get_id(),
                v->get_date() 
            );

    }
    
    // Update virus accounting
    p->n_viruses++;
    size_t n_viruses = p->n_viruses;

    if (n_viruses <= p->viruses.size())
        p->viruses[n_viruses - 1] = std::make_shared< Virus<TSeq> >(*v);
    else
        p->viruses.push_back(std::make_shared< Virus<TSeq> >(*v));

    n_viruses--;

    // Notice that both agent and date can be changed in this case
    // as only the sequence is a shared_ptr itself.
    #ifdef EPI_DEBUG
    if (n_viruses >= p->viruses.size())
    {
        throw std::logic_error(
            "[epi-debug]::default_add_virus Index for new virus out of range."
            );
    }
    #endif
    p->viruses[n_viruses]->set_agent(p, n_viruses);
    p->viruses[n_viruses]->set_date(m->today());

    #ifdef EPI_DEBUG
    m->get_db().today_variant.at(v->get_id()).at(p->state)++;
    #else
    m->get_db().today_variant[v->get_id()][p->state]++;
    #endif

}

template<typename TSeq>
inline void default_add_tool(Action<TSeq> & a, Model<TSeq> * m)
{

    Agent<TSeq> * p = a.agent;
    ToolPtr<TSeq> t = a.tool;

    CHECK_COALESCE_(a.new_state, t->status_init, p->get_state())
    CHECK_COALESCE_(a.queue, t->queue_init, QueueValues::NoOne)
    
    // Update tool accounting
    p->n_tools++;
    size_t n_tools = p->n_tools;

    if (n_tools <= p->tools.size())
        p->tools[n_tools - 1] = std::make_shared< Tool<TSeq> >(*t);
    else
        p->tools.push_back(std::make_shared< Tool<TSeq> >(*t));

    n_tools--;

    p->tools[n_tools]->set_date(m->today());
    p->tools[n_tools]->set_agent(p, n_tools);

    m->get_db().today_tool[t->get_id()][p->state]++;

}

template<typename TSeq>
inline void default_rm_virus(Action<TSeq> & a, Model<TSeq> * model)
{

    Agent<TSeq> * p    = a.agent;    
    VirusPtr<TSeq> & v = a.agent->viruses[a.virus->pos_in_agent];
    
    CHECK_COALESCE_(a.new_state, v->status_post, p->get_state())
    CHECK_COALESCE_(a.queue, v->queue_post, -QueueValues::Everyone)

    if (--p->n_viruses > 0)
    {
        // The new virus will change positions
        p->viruses[p->n_viruses]->pos_in_agent = v->pos_in_agent;
        std::swap(
            p->viruses[p->n_viruses],   // Moving to the end
            p->viruses[v->pos_in_agent] // Moving to the beginning
            );
    }
    
    // Calling the virus action over the removed virus
    v->post_recovery(model);

    return;

}

template<typename TSeq>
inline void default_rm_tool(Action<TSeq> & a, Model<TSeq> * /*m*/)
{

    Agent<TSeq> * p   = a.agent;    
    ToolPtr<TSeq> & t = a.agent->tools[a.tool->pos_in_agent];

    CHECK_COALESCE_(a.new_state, t->status_post, p->get_state())
    CHECK_COALESCE_(a.queue, t->queue_post, QueueValues::NoOne)

    if (--p->n_tools > 0)
    {
        p->tools[p->n_tools]->pos_in_agent = t->pos_in_agent;
        std::swap(
            p->tools[t->pos_in_agent],
            p->tools[p->n_tools]
            );
    }

    return;

}

template<typename TSeq>
inline void default_add_entity(Action<TSeq> & a, Model<TSeq> *)
{

    Agent<TSeq> *  p = a.agent;
    Entity<TSeq> * e = a.entity;

    CHECK_COALESCE_(a.new_state, e->status_post, p->get_state())
    CHECK_COALESCE_(a.queue, e->queue_post, QueueValues::NoOne)

    // Checking the agent and the entity are not linked
    if ((p->get_n_entities() > 0) && (e->size() > 0))
    {

        if (p->get_n_entities() > e->size()) // Slower search through the agent
        {
            for (size_t i = 0u; i < e->size(); ++i)
                if(e->operator[](i)->get_id() == p->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }
        else                                 // Slower search through the entity
        {
            for (size_t i = 0u; i < p->get_n_entities(); ++i)
                if(p->get_entity(i).get_id() == e->get_id())
                    throw std::logic_error("An entity cannot be reassigned to an agent.");
        }

        // It means that agent and entity were not associated.
    }

    // Adding the entity to the agent
    if (++p->n_entities <= p->entities.size())
    {

        p->entities[p->n_entities - 1]           = e->get_id();
        p->entities_locations[p->n_entities - 1] = e->n_agents;

    } else
    {
        p->entities.push_back(e->get_id());
        p->entities_locations.push_back(e->n_agents);
    }

    // Adding the agent to the entity
    // Adding the entity to the agent
    if (++e->n_agents <= e->agents.size())
    {

        e->agents[e->n_agents - 1]          = p->get_id();
        // Adjusted by '-1' since the list of entities in the agent just grew.
        e->agents_location[e->n_agents - 1] = p->n_entities - 1;

    } else
    {
        e->agents.push_back(p->get_id());
        e->agents_location.push_back(p->n_entities - 1);
    }

    // Today was the last modification
    // e->date_last_add_or_remove = m->today();
    
}

template<typename TSeq>
inline void default_rm_entity(Action<TSeq> & a, Model<TSeq> * m)
{
    
    Agent<TSeq> *  p = a.agent;    
    Entity<TSeq> * e = a.entity;
    size_t idx_agent_in_entity = a.idx_agent;
    size_t idx_entity_in_agent = a.idx_object;

    CHECK_COALESCE_(a.new_state, e->status_post, p->get_state())
    CHECK_COALESCE_(a.queue, e->queue_post, QueueValues::NoOne)

    if (--p->n_entities > 0)
    {

        // When we move the end entity to the new location, the 
        // moved entity needs to reflect the change, i.e., where the
        // entity will now be located in the agent
        size_t agent_location_in_last_entity  = p->entities_locations[p->n_entities];
        Entity<TSeq> * last_entity = &m->get_entities()[p->entities[p->n_entities]]; ///< Last entity of the agent

        // The end entity will be located where the removed was
        last_entity->agents_location[agent_location_in_last_entity] = idx_entity_in_agent;

        // We now make the swap
        std::swap(
            p->entities[p->n_entities],
            p->entities[idx_entity_in_agent]
        );

    }

    if (--e->n_agents > 0)
    {

        // When we move the end agent to the new location, the 
        // moved agent needs to reflect the change, i.e., where the
        // agent will now be located in the entity
        size_t entity_location_in_last_agent = e->agents_location[e->n_agents];
        Agent<TSeq> * last_agent  = &m->get_agents()[e->agents[e->n_agents]]; ///< Last agent of the entity

        // The end entity will be located where the removed was
        last_agent->entities_locations[entity_location_in_last_agent] = idx_agent_in_entity;

        // We now make the swap
        std::swap(
            e->agents[e->n_agents],
            e->agents[idx_agent_in_entity]
        );

    }

    // Setting the date of the last removal
    // e->date_last_add_or_remove = m->today();

    return;

}
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//agent-actions-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



template<typename TSeq>
inline Agent<TSeq>::Agent() {}

template<typename TSeq>
inline Agent<TSeq>::Agent(Agent<TSeq> && p) :
    model(p.model),
    neighbors(std::move(p.neighbors)),
    neighbors_locations(std::move(p.neighbors_locations)),
    n_neighbors(p.n_neighbors),
    entities(std::move(p.entities)),
    entities_locations(std::move(p.entities_locations)),
    n_entities(p.n_entities),
    state(p.state),
    status_prev(p.status_prev), 
    status_last_changed(p.status_last_changed),
    id(p.id),
    viruses(std::move(p.viruses)),  /// Needs to be adjusted
    n_viruses(p.n_viruses),
    tools(std::move(p.tools)), /// Needs to be adjusted
    n_tools(p.n_tools),
    add_virus_(std::move(p.add_virus_)),
    add_tool_(std::move(p.add_tool_)),
    add_entity_(std::move(p.add_entity_)),
    rm_virus_(std::move(p.rm_virus_)),
    rm_tool_(std::move(p.rm_tool_)),
    rm_entity_(std::move(p.rm_entity_)),
    action_counter(p.action_counter)
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus

    int loc = 0;
    for (auto & v : viruses)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        v->agent = this;
        v->pos_in_agent = loc++;

    }

    loc = 0;
    for (auto & t : tools)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        t->agent     = this;
        t->pos_in_agent = loc++;

    }
    
}

template<typename TSeq>
inline Agent<TSeq>::Agent(const Agent<TSeq> & p) :
    model(p.model),
    neighbors(p.neighbors),
    neighbors_locations(p.neighbors_locations),
    n_neighbors(p.n_neighbors),
    entities(p.entities),
    entities_locations(p.entities_locations),
    n_entities(p.n_entities),
    sampled_agents(0u),
    sampled_agents_n(0u),
    sampled_agents_left_n(0u),
    date_last_build_sample(-99)
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus
    viruses.resize(p.get_n_viruses(), nullptr);
    n_viruses = viruses.size();
    for (size_t i = 0u; i < n_viruses; ++i)
    {

        // Will create a copy of the virus, with the exeption of
        // the virus code
        viruses[i] = std::make_shared<Virus<TSeq>>(*p.viruses[i]);
        viruses[i]->set_agent(this, i);

    }

    tools.resize(p.get_n_tools(), nullptr);
    n_tools = tools.size();
    for (size_t i = 0u; i < n_tools; ++i)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        tools[i] = std::make_shared<Tool<TSeq>>(*p.tools[i]);
        tools[i]->set_agent(this, i);

    }

    add_virus_ = p.add_virus_;
    add_tool_  = p.add_tool_;
    rm_virus_  = p.rm_virus_;
    rm_tool_   = p.rm_tool_;
    
}

template<typename TSeq>
inline Agent<TSeq> & Agent<TSeq>::operator=(
    const Agent<TSeq> & other_agent
) 
{

    model = other_agent.model;

    neighbors = other_agent.neighbors;
    neighbors_locations = other_agent.neighbors_locations;
    n_neighbors = other_agent.n_neighbors;

    entities = other_agent.entities;
    entities_locations = other_agent.entities_locations;
    n_entities = other_agent.n_entities;

    sampled_agents.clear();
    sampled_agents_n = 0;
    sampled_agents_left_n = 0;
    date_last_build_sample = -99;

    // neighbors           = other_agent.neighbors;
    // entities            = other_agent.entities;
    // entities_locations  = other_agent.entities_locations;
    // n_entities          = other_agent.n_entities;
    state              = other_agent.state;
    status_prev         = other_agent.status_prev;
    status_last_changed = other_agent.status_last_changed;
    id                  = other_agent.id;
    
    // viruses             = other_agent.viruses;
    n_viruses           = other_agent.n_viruses;
    viruses.resize(n_viruses, nullptr);
    for (size_t i = 0u; i < n_viruses; ++i)
    {
        viruses[i] = std::make_shared<Virus<TSeq>>(*other_agent.viruses[i]);
        viruses[i]->set_agent(this, i);
    }
    
    // tools               = other_agent.tools;
    n_tools             = other_agent.n_tools;
    for (size_t i = 0u; i < n_tools; ++i)
    {
        tools[i] = std::make_shared<Tool<TSeq>>(*other_agent.tools[i]);
        tools[i]->set_agent(this, i);
    }

    add_virus_          = other_agent.add_virus_;
    add_tool_           = other_agent.add_tool_;
    add_entity_         = other_agent.add_entity_;
    rm_virus_           = other_agent.rm_virus_;
    rm_tool_            = other_agent.rm_tool_;
    rm_entity_          = other_agent.rm_entity_;
    action_counter      = other_agent.action_counter;
    
    return *this;
    
}

template<typename TSeq>
inline void Agent<TSeq>::add_tool(
    ToolPtr<TSeq> tool,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
) {

    // Checking the virus exists
    if (tool->get_id() >= static_cast<int>(model->get_db().get_n_tools()))
        throw std::range_error("The tool with id: " + std::to_string(tool->get_id()) + 
            " has not been registered. There are only " + std::to_string(model->get_n_tools()) + 
            " included in the model.");
    

    model->actions_add(
        this, nullptr, tool, nullptr, status_new, queue, add_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::add_tool(
    Tool<TSeq> tool,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{
    ToolPtr<TSeq> tool_ptr = std::make_shared< Tool<TSeq> >(tool);
    add_tool(tool_ptr, model, status_new, queue);
}

template<typename TSeq>
inline void Agent<TSeq>::add_virus(
    VirusPtr<TSeq> virus,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    // Checking the virus exists
    if (virus->get_id() >= static_cast<int>(model->get_db().get_n_variants()))
        throw std::range_error("The virus with id: " + std::to_string(virus->get_id()) + 
            " has not been registered. There are only " + std::to_string(model->get_n_variants()) + 
            " included in the model.");

    model->actions_add(
        this, virus, nullptr, nullptr, status_new, queue, add_virus_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::add_virus(
    Virus<TSeq> virus,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{
    VirusPtr<TSeq> virus_ptr = std::make_shared< Virus<TSeq> >(virus);
    add_virus(virus_ptr, model, status_new, queue);
}

template<typename TSeq>
inline void Agent<TSeq>::add_entity(
    Entity<TSeq> & entity,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (model != nullptr)
    {

        model->actions_add(
            this, nullptr, nullptr, &entity, status_new, queue, add_entity_, -1, -1
        );

    }
    else // If no model is passed, then we assume that we only need to add the
         // model entity
    {

        Action<TSeq> a(
                this, nullptr, nullptr, &entity, status_new, queue, add_entity_,
                -1, -1
            );

        default_add_entity(a, model); /* passing model makes nothing */

    }

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    epiworld_fast_uint tool_idx,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (tool_idx >= n_tools)
        throw std::range_error(
            "The Tool you want to remove is out of range. This Agent only has " +
            std::to_string(n_tools) + " tools."
        );

    model->actions_add(
        this, nullptr, tools[tool_idx], nullptr, status_new, queue, rm_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_tool(
    ToolPtr<TSeq> & tool,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (tool->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model->actions_add(
        this, nullptr, tool, nullptr, status_new, queue, rm_tool_, -1, -1
        );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    epiworld_fast_uint virus_idx,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{
    if (virus_idx >= n_viruses)
        throw std::range_error(
            "The Virus you want to remove is out of range. This Agent only has " +
            std::to_string(n_viruses) + " viruses."
        );
    else if (n_viruses == 0u)
        throw std::logic_error(
            "There is no virus to remove here!"
        );

    #ifdef EPI_DEBUG
    if (viruses[virus_idx]->pos_in_agent >= static_cast<int>(n_viruses))
    {
        throw std::logic_error(
            "[epi-debug]::rm_virus the position in the agent is wrong."
            );
    }
    #endif

    model->actions_add(
        this, viruses[virus_idx], nullptr, nullptr, status_new, queue,
        default_rm_virus<TSeq>, -1, -1
        );
    
}

template<typename TSeq>
inline void Agent<TSeq>::rm_virus(
    VirusPtr<TSeq> & virus,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (virus->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model->actions_add(
        this, virus, nullptr, nullptr, status_new, queue,
        default_rm_virus<TSeq>, -1, -1
        );


}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    epiworld_fast_uint entity_idx,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (entity_idx >= n_entities)
        throw std::range_error(
            "The Entity you want to remove is out of range. This Agent only has " +
            std::to_string(n_entities) + " entitites."
        );
    else if (n_entities == 0u)
        throw std::logic_error(
            "There is entity to remove here!"
        );

    model->actions_add(
        this, nullptr, nullptr, model->entities[entity_idx], status_new, queue, 
        default_rm_entity, entities_locations[entity_idx], entity_idx
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_entity(
    Entity<TSeq> & entity,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    // Looking for entity location in the agent
    int entity_idx = -1;
    for (size_t i = 0u; i < n_entities; ++i)
    {
        if (entities[i] == entity->get_id())
            entity_idx = i;
    }

    if (entity_idx == -1)
        throw std::logic_error(
            "The agent " + std::to_string(id) + " is not associated with entity \"" +
            entity.get_name() + "\"."
            );


    model->actions_add(
        this, nullptr, nullptr, entities[entity_idx], status_new, queue, 
        default_rm_entity, entities_locations[entity_idx], entity_idx
    );
}

template<typename TSeq>
inline void Agent<TSeq>::rm_agent_by_virus(
    epiworld_fast_uint virus_idx,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (status_new == -99)
        status_new = state;

    if (virus_idx >= n_viruses)
        throw std::range_error(
            std::string("The virus trying to remove the agent is out of range. ") +
            std::string("This agent has only ") + std::to_string(n_viruses) + 
            std::string(" and you are trying to remove virus # ") +
            std::to_string(virus_idx) + std::string(".")
            );

    // Removing viruses
    for (size_t i = 0u; i < n_viruses; ++i)
    {
        if (i != virus_idx)
            rm_virus(i, model);
    }

    // Changing state to new_state
    VirusPtr<TSeq> & v = viruses[virus_idx];
    epiworld_fast_int dead_state, dead_queue;
    v->get_state(nullptr, nullptr, &dead_state);
    v->get_queue(nullptr, nullptr, &dead_queue);

    if (queue != -99)
        dead_queue = queue;

    change_state(
        model,
        // Either preserve the current state or apply a new one
        (dead_state < 0) ? state : static_cast<epiworld_fast_uint>(dead_state),

        // By default, it will be removed from the queue... unless the user
        // says the contrary!
        (dead_queue == -99) ? QueueValues::NoOne : dead_queue
    );

}

template<typename TSeq>
inline void Agent<TSeq>::rm_agent_by_virus(
    VirusPtr<TSeq> & virus,
    Model<TSeq> * model,
    epiworld_fast_int status_new,
    epiworld_fast_int queue
)
{

    if (virus->get_agent() == nullptr)
        throw std::logic_error("The virus trying to remove the agent has no host.");

    if (virus->get_agent()->id != id)
        throw std::logic_error("Viruses can only remove their hosts'.");

    rm_agent_by_virus(
        virus->pos_in_agent,
        model,
        status_new,
        queue
    );

}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {

    return model->susceptibility_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_transmission_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->transmission_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_recovery_enhancer(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->recovery_enhancer_mixer(this, v, model);
}

template<typename TSeq>
inline epiworld_double Agent<TSeq>::get_death_reduction(
    VirusPtr<TSeq> v,
    Model<TSeq> * model
) {
    return model->death_reduction_mixer(this, v, model);
}

template<typename TSeq>
inline int Agent<TSeq>::get_id() const
{
    return id;
}

template<typename TSeq>
inline Viruses<TSeq> Agent<TSeq>::get_viruses() {

    return Viruses<TSeq>(*this);

}

template<typename TSeq>
inline const Viruses_const<TSeq> Agent<TSeq>::get_viruses() const {

    return Viruses_const<TSeq>(*this);
    
}

template<typename TSeq>
inline VirusPtr<TSeq> & Agent<TSeq>::get_virus(int i) {
    return viruses.at(i);
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_viruses() const noexcept
{
    return n_viruses;
}

template<typename TSeq>
inline Tools<TSeq> Agent<TSeq>::get_tools() {
    return Tools<TSeq>(*this);
}

template<typename TSeq>
inline const Tools_const<TSeq> Agent<TSeq>::get_tools() const {
    return Tools_const<TSeq>(*this);
}

template<typename TSeq>
inline ToolPtr<TSeq> & Agent<TSeq>::get_tool(int i)
{
    return tools.at(i);
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_tools() const noexcept
{
    return n_tools;
}

template<typename TSeq>
inline void Agent<TSeq>::mutate_variant()
{

    for (auto & v : viruses)
        v->mutate();

}

template<typename TSeq>
inline void Agent<TSeq>::add_neighbor(
    Agent<TSeq> & p,
    bool check_source,
    bool check_target
) {
    // Can we find the neighbor?
    bool found = false;
    if (check_source)
    {

        for (auto & n: neighbors)    
            if (static_cast<int>(n) == p.get_id())
            {
                found = true;
                break;
            }

    }

    // Three things going on here:
    // - Where in the neighbor will this be
    // - What is the neighbor's id
    // - Increasing the number of neighbors
    if (!found)
    {

        neighbors_locations.push_back(p.get_n_neighbors());
        neighbors.push_back(p.get_id());
        n_neighbors++;

    }


    found = false;
    if (check_target)
    {

        for (auto & n: p.neighbors)
            if (static_cast<int>(n) == id)
            {
                found = true;
                break;
            }
    
    }

    if (!found)
    {

        p.neighbors_locations.push_back(n_neighbors - 1);
        p.neighbors.push_back(id);
        p.n_neighbors++;
        
    }
    

}

template<typename TSeq>
inline void Agent<TSeq>::swap_neighbors(
    Agent<TSeq> & other,
    size_t n_this,
    size_t n_other
)
{

    // Getting the agents
    auto & pop = model->population;
    auto & neigh_this  = pop[neighbors[n_this]];
    auto & neigh_other = pop[other.neighbors[n_other]];

    // Getting the locations in the neighbors
    size_t loc_this_in_neigh = neighbors_locations[n_this];
    size_t loc_other_in_neigh = other.neighbors_locations[n_other];

    // Changing ids
    std::swap(neighbors[n_this], other.neighbors[n_other]);

    if (!model->directed)
    {
        std::swap(
            neigh_this.neighbors[loc_this_in_neigh],
            neigh_other.neighbors[loc_other_in_neigh]
            );

        // Changing the locations
        std::swap(neighbors_locations[n_this], other.neighbors_locations[n_other]);
        
        std::swap(
            neigh_this.neighbors_locations[loc_this_in_neigh],
            neigh_other.neighbors_locations[loc_other_in_neigh]
            );
    }

}

template<typename TSeq>
inline std::vector< Agent<TSeq> *> Agent<TSeq>::get_neighbors()
{
    std::vector< Agent<TSeq> * > res(n_neighbors, nullptr);
    for (size_t i = 0u; i < n_neighbors; ++i)
        res[i] = &model->population[neighbors[i]];

    return res;
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_neighbors() const
{
    return n_neighbors;
}

template<typename TSeq>
inline void Agent<TSeq>::change_state(
    Model<TSeq> * model,
    epiworld_fast_uint new_state,
    epiworld_fast_int queue
    )
{

    model->actions_add(
        this, nullptr, nullptr, nullptr, new_state, queue, nullptr, -1, -1
    );
    
    return;

}

template<typename TSeq>
inline const epiworld_fast_uint & Agent<TSeq>::get_state() const {
    return state;
}

template<typename TSeq>
inline void Agent<TSeq>::reset()
{

    this->viruses.clear();
    n_viruses = 0u;

    this->tools.clear();
    n_tools = 0u;

    this->state = 0u;
    this->status_prev = 0u;

    this->status_last_changed = -1;
    
}

template<typename TSeq>
inline bool Agent<TSeq>::has_tool(epiworld_fast_uint t) const
{

    for (auto & tool : tools)
        if (tool->get_id() == t)
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_tool(std::string name) const
{

    for (auto & tool : tools)
        if (tool->get_name() == name)
            return true;

    return false;

}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(epiworld_fast_uint t) const
{
    for (auto & v : viruses)
        if (v->get_id() == t)
            return true;

    return false;
}

template<typename TSeq>
inline bool Agent<TSeq>::has_virus(std::string name) const
{
    
    for (auto & v : viruses)
        if (v->get_name() == name)
            return true;

    return false;

}

template<typename TSeq>
inline void Agent<TSeq>::print(
    Model<TSeq> * model,
    bool compressed
    ) const
{

    if (compressed)
    {
        printf_epiworld(
            "Agent: %i, state: %s (%lu), Nvirus: %lu, NTools: %lu, NNeigh: %lu\n",
            id, model->states_labels[state].c_str(), state, n_viruses, n_tools, neighbors.size()
        );
    }
    else {

        printf_epiworld("Information about agent id %i\n", this->id);
        printf_epiworld("  State        : %s (%lu)\n", model->states_labels[state].c_str(), state);
        printf_epiworld("  Virus count  : %lu\n", n_viruses);
        printf_epiworld("  Tool count   : %lu\n", n_tools);
        printf_epiworld("  Neigh. count : %lu\n", neighbors.size());

        size_t nfeats = model->get_agents_data_ncols();
        if (nfeats > 0)
        {

            printf_epiworld("This model includes features (%lu): [ ", nfeats);

            int max_to_show = static_cast<int>((nfeats > 10)? 10 : nfeats);

            for (int k = 0; k < max_to_show; ++k)
            {
                printf_epiworld("%.2f", this->operator[](k));

                if (k != (max_to_show - 1))
                {
                    printf_epiworld(", ");
                } else {
                    printf_epiworld(" ]\n");
                }

            }
            
        }

    }

    return;

}

template<typename TSeq>
inline double & Agent<TSeq>::operator()(size_t j)
{

    if (model->agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model->agents_data + j * model->size() + id);

}

template<typename TSeq>
inline double & Agent<TSeq>::operator[](size_t j)
{
    return *(model->agents_data + j * model->size() + id);
}

template<typename TSeq>
inline double Agent<TSeq>::operator()(size_t j) const
{

    if (model->agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model->agents_data + j * model->size() + id);

}

template<typename TSeq>
inline double Agent<TSeq>::operator[](size_t j) const
{
    return *(model->agents_data + j * model->size() + id);
}

template<typename TSeq>
inline Entities<TSeq> Agent<TSeq>::get_entities()
{
    return Entities<TSeq>(*this);
}

template<typename TSeq>
inline const Entities_const<TSeq> Agent<TSeq>::get_entities() const
{
    return Entities_const<TSeq>(*this);
}

template<typename TSeq>
inline const Entity<TSeq> & Agent<TSeq>::get_entity(size_t i) const
{
    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->entities[entities[i]];
}

template<typename TSeq>
inline Entity<TSeq> & Agent<TSeq>::get_entity(size_t i)
{
    if (i >= n_entities)
        throw std::range_error("Trying to get to an agent's entity outside of the range.");

    return model->entities[entities[i]];
}

template<typename TSeq>
inline size_t Agent<TSeq>::get_n_entities() const
{
    return n_entities;
}

template<typename TSeq>
inline bool Agent<TSeq>::operator==(const Agent<TSeq> & other) const
{

    EPI_DEBUG_FAIL_AT_TRUE(
        n_neighbors != other.n_neighbors,
        "Agent:: n_eighbors don't match"
        )

    
    for (size_t i = 0u; i < n_neighbors; ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            neighbors[i] != other.neighbors[i],
            "Agent:: neighbor[i] don't match"
        )
    }
    
    EPI_DEBUG_FAIL_AT_TRUE(
        n_entities != other.n_entities,
        "Agent:: n_entities don't match"
        )
    
    
    for (size_t i = 0u; i < n_entities; ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            entities[i] != other.entities[i],
            "Agent:: entities[i] don't match"
        )
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        state != other.state,
        "Agent:: state don't match"
        )
        

    EPI_DEBUG_FAIL_AT_TRUE(
        status_prev != other.status_prev,
        "Agent:: status_prev don't match"
        )
        

    // EPI_DEBUG_FAIL_AT_TRUE(
    //     status_last_changed != other.status_last_changed,
    //     "Agent:: status_last_changed don't match"
    //     ) ///< Last time the agent was updated.
        
    
    EPI_DEBUG_FAIL_AT_TRUE(
        n_viruses != other.n_viruses,
        "Agent:: n_viruses don't match"
        )
        

    for (size_t i = 0u; i < n_viruses; ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            *viruses[i] != *other.viruses[i],
            "Agent:: viruses[i] don't match"
        )
         
    }

    EPI_DEBUG_FAIL_AT_TRUE(n_tools != other.n_tools, "Agent:: n_tools don't match")

    for (size_t i = 0u; i < n_tools; ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            tools[i] != other.tools[i],
            "Agent:: tools[i] don't match"
        )
         
    }   
    
    return true;
    
}

#undef CHECK_COALESCE_

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/agent-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/agentssample-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_AGENTS_BONES_HPP
#define EPIWORLD_AGENTS_BONES_HPP

class SAMPLETYPE {
public:
    static const int MODEL  = 0;
    static const int ENTITY = 1;
    static const int AGENT  = 2;
};

// template<typename TSeq>
// class Agent;

// template<typename TSeq>
// class Model;

// template<typename TSeq>
// class Entity;

/**
 * @brief Sample of agents
 * 
 * This class allows sampling agents from Entity<TSeq> and Model<TSeq>.
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class AgentsSample {
private:

    size_t sample_size = 0u;

    std::vector< Agent<TSeq>* > * agents = nullptr; ///< Pointer to sample of agents
    size_t * agents_n = nullptr;                    ///< Size of sample of agents
    
    std::vector< size_t > * agents_left = nullptr;  ///< Pointer to agents left (iota)
    size_t * agents_left_n = nullptr;               ///< Size of agents left

    Model<TSeq> * model   = nullptr;   ///< Extracts runif() and (if the case) population.
    Entity<TSeq> * entity = nullptr; ///
    Agent<TSeq> * agent   = nullptr;
    
    int sample_type = SAMPLETYPE::AGENT;

    void sample_n(size_t n); ///< Backbone function for sampling


public:

    // Not available (for now)
    AgentsSample() = delete;                       ///< Default constructor
    AgentsSample(const AgentsSample<TSeq> & a) = delete; ///< Copy constructor
    AgentsSample(AgentsSample<TSeq> && a) = delete;      ///< Move constructor

    AgentsSample(Model<TSeq> & model_, size_t n, bool truncate = false);
    AgentsSample(Model<TSeq> * model, Entity<TSeq> & entity_, size_t n, bool truncate = false);
    AgentsSample(Model<TSeq> * model, Agent<TSeq> & agent_, size_t n, bool truncate = false);

    ~AgentsSample();

    typename std::vector< Agent<TSeq> * >::iterator begin();
    typename std::vector< Agent<TSeq> * >::iterator end();

    Agent<TSeq> * operator[](size_t n);
    Agent<TSeq> * operator()(size_t n);
    size_t size() const noexcept;

};

template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> & model_,
    size_t n,
    bool truncate
    ) {

    if (truncate)
    {
        
        if (n > model_.size())
            n = model_.size();

    } else if (n > model_.size())
        throw std::logic_error(
            "There are only " + std::to_string(model_.size()) + " agents. You cannot " +
            "sample " + std::to_string(n));

    sample_size = n;
    model       = &model_;
    sample_type = SAMPLETYPE::MODEL;

    agents   = &model_.sampled_population;
    agents_n = &model_.sampled_population_n;

    agents_left   = &model_.population_left;
    agents_left_n = &model_.population_left_n;

    sample_n(n);
    
    return; 

}

template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> * model,
    Entity<TSeq> & entity_,
    size_t n, bool truncate
    ) {

    if (truncate)
    {

        if (n > entity_.size())
            n = entity_.size();

    } else if (n > entity_.size())
        throw std::logic_error(
            "There are only " + std::to_string(entity_.size()) + " agents. You cannot " +
            "sample " + std::to_string(n));

    sample_size = n;
    model       = &entity_.model;
    sample_type = SAMPLETYPE::ENTITY;

    agents   = &entity_.sampled_agents;
    agents_n = &entity_.sampled_agents_n;

    agents_left   = &entity_.sampled_agents_left;
    agents_left_n = &entity_.sampled_agents_left_n;

    sample_n(n);

    return; 

}

/**
 * @brief Sample from the agent's entities
 * 
 * For example, how many individuals the agent contacts in a given point in time.
 * 
 * @tparam TSeq 
 * @param agent_ 
 * @param n Sample size
 * @param truncate If the agent has fewer than `n` connections, then truncate = true
 * will automatically reduce the number of possible samples. Otherwise, if false, then
 * it returns an error.
 */
template<typename TSeq>
inline AgentsSample<TSeq>::AgentsSample(
    Model<TSeq> * model,
    Agent<TSeq> & agent_,
    size_t n,
    bool truncate
    )
{

    sample_type = SAMPLETYPE::AGENT;
    
    agent = &agent_;

    agents   = &agent_.sampled_agents;
    agents_n = &agent_.sampled_agents_n;

    agents_left   = &agent_.sampled_agents_left;
    agents_left_n = &agent_.sampled_agents_left_n;

    // Computing the cumulative sum of counts across entities
    size_t agents_in_entities = 0;
    Entities<TSeq> entities_a = agent->get_entities();

    std::vector< int > cum_agents_count(entities_a.size(), 0);
    int idx = -1;
    for (auto & e : entities_a)
    {
        if (++idx == 0)
            cum_agents_count[idx] = (e->size() - 1u);
        else
            cum_agents_count[idx] = (
                (e->size() - 1u) + 
                cum_agents_count[idx - 1]
            );

        agents_in_entities += (e->size() - 1u);
    }

    if (truncate)
    {
        
        if (n > agents_in_entities)
            n = agents_in_entities;

    } else if (n > agents_in_entities)
        throw std::logic_error(
            "There are only " + std::to_string(agents_in_entities) + " agents. You cannot " +
            "sample " + std::to_string(n));

    sample_size = n;

    if (agents->size() < n)
        agents->resize(n);

    for (size_t i = 0u; i < n; ++i)
    {
        int jth = std::floor(model->runif() * agents_in_entities);
        for (size_t e = 0u; e < cum_agents_count.size(); ++e)
        {
            // Are we in the limit?
            if (jth <= cum_agents_count[e])
            {
                if (e == 0) // From the first group
                    agents->operator[](i) = entities_a[e][jth];
                else
                    agents->operator[](i) = entities_a[e][jth - cum_agents_count[e - 1]];
                
                break;
            }

        }
    }

    return; 

}

template<typename TSeq>
inline AgentsSample<TSeq>::~AgentsSample() {}

template<typename TSeq>
inline size_t AgentsSample<TSeq>::size() const noexcept
{ 
    return this->sample_size;
}

template<typename TSeq>
inline Agent<TSeq> * AgentsSample<TSeq>::operator[](size_t i)
{

    return agents->operator[](i);

}

template<typename TSeq>
inline Agent<TSeq> * AgentsSample<TSeq>::operator()(size_t i)
{

    if (i >= this->sample_size)
        throw std::range_error("The requested agent is out of range.");

    return agents->operator[](i);

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator AgentsSample<TSeq>::begin()
{

    if (sample_size > 0u)
        return agents->begin();
    else
        return agents->end();

}

template<typename TSeq>
inline typename std::vector< Agent<TSeq> * >::iterator AgentsSample<TSeq>::end()
{

    return agents->begin() + sample_size;

}

template<typename TSeq>
inline void AgentsSample<TSeq>::sample_n(size_t n)
{

    // Checking if the size of the entity has changed (or hasn't been initialized)
    if (sample_type == SAMPLETYPE::MODEL)
    {

        if (model->size() != agents_left->size())
        {
            agents_left->resize(model->size(), 0u);
            std::iota(agents_left->begin(), agents_left->end(), 0u);
        }


    } else if (sample_type == SAMPLETYPE::ENTITY) {

        if (entity->size() != agents_left->size())
        {

            agents_left->resize(entity->size(), 0u);
            std::iota(agents_left->begin(), agents_left->end(), 0u);

        }

    } 

    // Restart the counter of agents left
    *agents_left_n = agents_left->size();

    if (agents->size() < sample_size)
        agents->resize(sample_size, nullptr);

    if (sample_type == SAMPLETYPE::MODEL)
    {

        for (size_t i = 0u; i < n; ++i)
        {

            size_t ith = agents_left->operator[](model->runif() * ((*agents_left_n)--));
            agents->operator[](i) = &model->population[ith];

            // Updating list
            std::swap(agents_left->operator[](ith), agents_left->operator[](*agents_left_n));

        }

    } else if (sample_type == SAMPLETYPE::ENTITY) {

        for (size_t i = 0u; i < n; ++i)
        {

            size_t ith = agents_left->operator[](model->runif() * (--(*agents_left_n)));
            agents->operator[](i) = entity->agents[ith];

            // Updating list
            std::swap(agents_left->operator[](ith), agents_left->operator[](*agents_left_n));

        }

    }

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/agentssample-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld/models/models.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_HPP
#define EPIWORLD_MODELS_HPP

namespace epimodels {

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/sis.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SIS_HPP 
#define EPIWORLD_MODELS_SIS_HPP

/**
 * @brief Template for a Susceptible-Infected-Susceptible (SIS) model
 * 
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery rate of the immune system
 */
template<typename TSeq = int>
class ModelSIS : public epiworld::Model<TSeq>
{

public:

    ModelSIS() {};

    ModelSIS(
        ModelSIS<TSeq> & model,
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double recovery
    );

    ModelSIS(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double recovery
    );

};

template<typename TSeq>
inline ModelSIS<TSeq>::ModelSIS(
    ModelSIS<TSeq> & model,
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double recovery
    )
{

    model.set_name("Susceptible-Infected-Susceptible (SIS)");

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Infected", epiworld::default_update_exposed<TSeq>);

    // Setting up parameters
    model.add_param(infectiousness, "Infection rate");
    model.add_param(recovery, "Recovery rate");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1,0,0);
    
    virus.set_prob_infecting(&model("Infection rate"));
    virus.set_prob_recovery(&model("Recovery rate"));
    virus.set_prob_death(0.0);
    
    model.add_virus(virus, prevalence);

    return;

}

template<typename TSeq>
inline ModelSIS<TSeq>::ModelSIS(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double recovery
    )
{

    ModelSIS<TSeq>(
        *this,
        vname,
        prevalence,
        infectiousness,
        recovery
    );    

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/sis.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/sir.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_SIR_H 
#define EPIWORLD_SIR_H

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery rate of the immune system
 */
template<typename TSeq = int>
class ModelSIR : public epiworld::Model<TSeq>
{
public:

    ModelSIR() {};

    ModelSIR(
        ModelSIR<TSeq> & model,
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double recovery
    );

    ModelSIR(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double recovery
    );
    
};

template<typename TSeq>
inline ModelSIR<TSeq>::ModelSIR(
    ModelSIR<TSeq> & model,
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double recovery
    )
{

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Infected", epiworld::default_update_exposed<TSeq>);
    model.add_state("Recovered");

    // Setting up parameters
    model.add_param(recovery, "Prob. of Recovery");
    model.add_param(infectiousness, "Infectiousness");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1,2,2);
    
    virus.set_prob_recovery(&model("Prob. of Recovery"));
    virus.set_prob_infecting(&model("Infectiousness"));
    
    model.add_virus(virus, prevalence);

    model.set_name("Susceptible-Infected-Recovered (SIR)");

    return;
   
}

template<typename TSeq>
inline ModelSIR<TSeq>::ModelSIR(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double recovery
    )
{

    ModelSIR<TSeq>(
        *this,
        vname,
        prevalence,
        infectiousness,
        recovery
        );

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/sir.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/seir.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SEIR_HPP
#define EPIWORLD_MODELS_SEIR_HPP

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery rate of the immune system
 */
template<typename TSeq = int>
class ModelSEIR : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int REMOVED     = 3;

public:

    ModelSEIR() {};

    ModelSEIR(
        ModelSEIR<TSeq> & model,
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double incubation_days,
        epiworld_double recovery
    );

    ModelSEIR(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double infectiousness,
        epiworld_double incubation_days,
        epiworld_double recovery
    );
    
    epiworld::UpdateFun<TSeq> update_exposed_seir = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Model<TSeq> * m
    ) -> void {
        // Does the agent become infected?
        if (m->runif() < 1.0/(m->par("Incubation days")))
            p->change_state(m, ModelSEIR<TSeq>::INFECTED);

        return;    
    };
      

    epiworld::UpdateFun<TSeq> update_infected_seir = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Model<TSeq> * m
    ) -> void {
        // Does the agent recover?
        if (m->runif() < (m->par("Immune recovery")))
            p->rm_virus(0, m);

        return;    
    };

};


template<typename TSeq>
inline ModelSEIR<TSeq>::ModelSEIR(
    ModelSEIR<TSeq> & model,
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double incubation_days,
    epiworld_double recovery
    )
{

    // Adding statuses
    model.add_state("Susceptible", epiworld::default_update_susceptible<TSeq>);
    model.add_state("Exposed", model.update_exposed_seir);
    model.add_state("Infected", model.update_infected_seir);
    model.add_state("Removed");

    // Setting up parameters
    model.add_param(infectiousness, "Infectiousness");
    model.add_param(incubation_days, "Incubation days");
    model.add_param(recovery, "Immune recovery");

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(ModelSEIR<TSeq>::EXPOSED, ModelSEIR<TSeq>::REMOVED, ModelSEIR<TSeq>::REMOVED);

    virus.set_prob_infecting(&model("Infectiousness"));
    
    // Adding the tool and the virus
    model.add_virus(virus, prevalence);
    
    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR)");

    return;
   
}

template<typename TSeq>
inline ModelSEIR<TSeq>::ModelSEIR(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double infectiousness,
    epiworld_double incubation_days,
    epiworld_double recovery
    )
{

    ModelSEIR<TSeq>(
        *this,
        vname,
        prevalence,
        infectiousness,
        incubation_days,
        recovery
        );

    return;

}



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/seir.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/surveillance.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SURVEILLANCE_HPP
#define EPIWORLD_MODELS_SURVEILLANCE_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSURV : public epiworld::Model<TSeq> {

private:
    // Status
    static const int SUSCEPTIBLE           = 0;
    static const int LATENT                = 1;
    static const int SYMPTOMATIC           = 2;
    static const int SYMPTOMATIC_ISOLATED  = 3; // sampled and discovered
    static const int ASYMPTOMATIC          = 4;
    static const int ASYMPTOMATIC_ISOLATED = 5;
    static const int RECOVERED             = 6;
    static const int REMOVED               = 7;

public:

    /**
     * @name Construct a new ModelSURV object
     * 
     * The ModelSURV class simulates a survaillence model where agents can be
     * isolated, even if asyptomatic.
     * 
     * @param vname String. Name of the virus
     * @param prevalence Integer. Number of initial cases of the virus.
     * @param efficacy_vax Double. Efficacy of the vaccine (1 - P(acquire the disease)).
     * @param latent_period Double. Shape parameter of a `Gamma(latent_period, 1)`
     *   distribution. This coincides with the expected number of latent days.
     * @param infect_period Double. Shape parameter of a `Gamma(infected_period, 1)`
     *   distribution. This coincides with the expected number of infectious days.
     * @param prob_symptoms Double. Probability of generating symptoms.
     * @param prop_vaccinated Double. Probability of vaccination. Coincides with
     *   the initial prevalence of vaccinated individuals.
     * @param prop_vax_redux_transm Double. Factor by which the vaccine reduces
     *   transmissibility.
     * @param prop_vax_redux_infect Double. Factor by which the vaccine reduces
     *   the chances of becoming infected.
     * @param surveillance_prob Double. Probability of testing an agent.
     * @param prob_transmission Double. Raw transmission probability.
     * @param prob_death Double. Raw probability of death for symptomatic individuals.
     * @param prob_noreinfect Double. Probability of no re-infection.
     * 
     * @details
     * This model features the following states:
     * 
     * - Susceptible
     * - Latent
     * - Symptomatic
     * - Symptomatic isolated
     * - Asymptomatic
     * - Asymptomatic isolated
     * - Recovered
     * - Removed    
     * 
     * @returns An object of class `epiworld_surv`
     * 
     */
    ///@{
    ModelSURV() {};

    ModelSURV(
        ModelSURV<TSeq> & model,
        std::string vname,
        epiworld_fast_uint prevalence               = 50,
        epiworld_double efficacy_vax          = 0.9,
        epiworld_double latent_period         = 3u,
        epiworld_double infect_period         = 6u,
        epiworld_double prob_symptoms         = 0.6,
        epiworld_double prop_vaccinated       = 0.25,
        epiworld_double prop_vax_redux_transm = 0.5,
        epiworld_double prop_vax_redux_infect = 0.5,
        epiworld_double surveillance_prob     = 0.001,
        epiworld_double prob_transmission     = 1.0,
        epiworld_double prob_death            = 0.001,
        epiworld_double prob_noreinfect       = 0.9
    );

    ModelSURV(
        std::string vname,
        epiworld_fast_uint prevalence         = 50,
        epiworld_double efficacy_vax          = 0.9,
        epiworld_double latent_period         = 3u,
        epiworld_double infect_period         = 6u,
        epiworld_double prob_symptoms         = 0.6,
        epiworld_double prop_vaccinated       = 0.25,
        epiworld_double prop_vax_redux_transm = 0.5,
        epiworld_double prop_vax_redux_infect = 0.5,
        epiworld_double surveillance_prob     = 0.001,
        epiworld_double prob_transmission     = 1.0,
        epiworld_double prob_death            = 0.001,
        epiworld_double prob_noreinfect       = 0.9
    );
    ///@}

};

template<typename TSeq>
inline ModelSURV<TSeq>::ModelSURV(
    ModelSURV<TSeq> & model,
    std::string vname,
    epiworld_fast_uint prevalence,
    epiworld_double efficacy_vax,
    epiworld_double latent_period,
    epiworld_double infect_period,
    epiworld_double prob_symptoms,
    epiworld_double prop_vaccinated,
    epiworld_double prop_vax_redux_transm,
    epiworld_double prop_vax_redux_infect,
    epiworld_double surveillance_prob,
    epiworld_double prob_transmission,
    epiworld_double prob_death,
    epiworld_double prob_noreinfect
    )
{

    EPI_NEW_UPDATEFUN_LAMBDA(surveillance_update_susceptible, TSeq) {

        // This computes the prob of getting any neighbor variant
        epiworld_fast_uint nvariants_tmp = 0u;
        for (auto & neighbor: p->get_neighbors()) 
        {
                    
            for (size_t i = 0u; i < neighbor->get_n_viruses(); ++i) 
            { 

                auto & v = neighbor->get_virus(i);
                    
                /* And it is a function of susceptibility_reduction as well */ 
                epiworld_double tmp_transmission = 
                    (1.0 - p->get_susceptibility_reduction(v, m)) * 
                    v->get_prob_infecting(m) * 
                    (1.0 - neighbor->get_transmission_reduction(v, m)) 
                    ; 
            
                m->array_double_tmp[nvariants_tmp]  = tmp_transmission;
                m->array_virus_tmp[nvariants_tmp++] = &(*v);
                
            } 
        }

        // No virus to compute on
        if (nvariants_tmp == 0)
            return;

        // Running the roulette
        int which = roulette(nvariants_tmp, m);

        if (which < 0)
            return;

        p->add_virus(*m->array_virus_tmp[which], m); 
        return;

    };


    epiworld::UpdateFun<TSeq> surveillance_update_exposed = 
    [](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
    {

        epiworld::VirusPtr<TSeq> & v = p->get_virus(0u); 
        epiworld_double p_die = v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, m)); 
        
        epiworld_fast_uint days_since_exposed = m->today() - v->get_date();
        epiworld_fast_uint status = p->get_state();

        // Figuring out latent period
        if (v->get_data().size() == 0u)
        {
            epiworld_double latent_days = m->rgamma(m->par("Latent period"), 1.0);
            v->get_data().push_back(latent_days);

            v->get_data().push_back(
                m->rgamma(m->par("Infect period"), 1.0) + latent_days
            );
        }
        
        // If still latent, nothing happens
        if (days_since_exposed <= v->get_data()[0u])
            return;

        // If past days infected + latent, then bye.
        if (days_since_exposed >= v->get_data()[1u])
        {
            p->rm_virus(0, m);
            return;
        }

        // If it is infected, then it can be asymptomatic or symptomatic
        if (status == ModelSURV<TSeq>::LATENT)
        {

            // Will be symptomatic?
            if (EPI_RUNIF() < m->par("Prob of symptoms"))
                p->change_state(m, ModelSURV<TSeq>::SYMPTOMATIC);
            else
                p->change_state(m, ModelSURV<TSeq>::ASYMPTOMATIC);
            
            return;

        }
        
        // Otherwise, it can be removed
        if (EPI_RUNIF() < p_die)
        {
            p->change_state(m, ModelSURV<TSeq>::REMOVED, -1);
            return;
        }
        
        return;

    };

    std::vector< epiworld_fast_uint > exposed_state = {
        SYMPTOMATIC,
        SYMPTOMATIC_ISOLATED,
        ASYMPTOMATIC,
        ASYMPTOMATIC_ISOLATED,
        LATENT
    };

    epiworld::GlobalFun<TSeq> surveillance_program = 
    [exposed_state](
        epiworld::Model<TSeq>* m
        ) -> void
    {

        // How many will we find
        std::binomial_distribution<> bdist(m->size(), m->par("Surveilance prob."));
        int nsampled = bdist(m->get_rand_endgine());

        int to_go = nsampled + 1;

        epiworld_double ndetected        = 0.0;
        epiworld_double ndetected_asympt = 0.0;
        
        auto & pop = m->get_agents();
        std::vector< bool > sampled(m->size(), false);
        
        while (to_go-- > 0)
        {

            // Who is the lucky one
            epiworld_fast_uint i = static_cast<epiworld_fast_uint>(std::floor(EPI_RUNIF() * m->size()));

            if (sampled[i])
                continue;

            sampled[i] = true;
            epiworld::Agent<TSeq> * p = &pop[i];
            
            // If still exposed for the next term
            if (epiworld::IN(p->get_state(), exposed_state ))
            {

                ndetected += 1.0;
                if (p->get_state() == ModelSURV<TSeq>::ASYMPTOMATIC)
                {
                    ndetected_asympt += 1.0;
                    p->change_state(m, ModelSURV<TSeq>::ASYMPTOMATIC_ISOLATED);
                }
                else 
                {
                    p->change_state(m, ModelSURV<TSeq>::SYMPTOMATIC_ISOLATED);
                }

            }

        }

        // Writing the user data
        std::vector< int > totals;
        m->get_db().get_today_total(nullptr,&totals);
        m->add_user_data(
            {
                static_cast<epiworld_double>(nsampled),
                ndetected,
                ndetected_asympt,
                static_cast<epiworld_double>(totals[ModelSURV<TSeq>::ASYMPTOMATIC])
            }
            );


    };

    model.add_state("Susceptible", surveillance_update_susceptible);
    model.add_state("Latent", surveillance_update_exposed);
    model.add_state("Symptomatic", surveillance_update_exposed);
    model.add_state("Symptomatic isolated", surveillance_update_exposed);
    model.add_state("Asymptomatic", surveillance_update_exposed);
    model.add_state("Asymptomatic isolated", surveillance_update_exposed);
    model.add_state("Recovered");
    model.add_state("Removed");

    // General model parameters
    model.add_param(latent_period, "Latent period");
    model.add_param(infect_period, "Infect period");
    model.add_param(prob_symptoms, "Prob of symptoms");
    model.add_param(surveillance_prob, "Surveilance prob.");
    model.add_param(efficacy_vax, "Vax efficacy");
    model.add_param(prop_vax_redux_transm, "Vax redux transmission");
    model.add_param(prob_transmission, "Prob of transmission");
    model.add_param(prob_death, "Prob. death");
    model.add_param(prob_noreinfect, "Prob. no reinfect");

    // Virus ------------------------------------------------------------------
    epiworld::Virus<TSeq> covid("Covid19");
    covid.set_state(LATENT, RECOVERED, REMOVED);
    covid.set_post_immunity(&model("Prob. no reinfect"));
    covid.set_prob_death(&model("Prob. death"));

    epiworld::VirusFun<TSeq> ptransmitfun = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Virus<TSeq> & v,
        epiworld::Model<TSeq> * m
        ) -> epiworld_double
    {
        // No chance of infecting
        epiworld_fast_uint  s = p->get_state();
        if (s == ModelSURV<TSeq>::LATENT)
            return static_cast<epiworld_double>(0.0);
        else if (s == ModelSURV<TSeq>::SYMPTOMATIC_ISOLATED)
            return static_cast<epiworld_double>(0.0);
        else if (s == ModelSURV<TSeq>::ASYMPTOMATIC_ISOLATED)
            return static_cast<epiworld_double>(0.0);

        // Otherwise
        return m->par("Prob of transmission");
    };

    covid.set_prob_infecting_fun(ptransmitfun);
    
    model.add_virus_n(covid, prevalence);

    model.set_user_data({"nsampled", "ndetected", "ndetected_asympt", "nasymptomatic"});
    model.add_global_action(surveillance_program,-1);
   
    // Vaccine tool -----------------------------------------------------------
    epiworld::Tool<TSeq> vax("Vaccine");
    vax.set_susceptibility_reduction(&model("Vax efficacy"));
    vax.set_transmission_reduction(&model("Vax redux transmission"));
    
    model.add_tool(vax, prop_vaccinated);

    model.set_name("Surveillance");

    return;

}

template<typename TSeq>
inline ModelSURV<TSeq>::ModelSURV(
    std::string vname,
    epiworld_fast_uint prevalence,
    epiworld_double efficacy_vax,
    epiworld_double latent_period,
    epiworld_double infect_period,
    epiworld_double prob_symptoms,
    epiworld_double prop_vaccinated,
    epiworld_double prop_vax_redux_transm,
    epiworld_double prop_vax_redux_infect,
    epiworld_double surveillance_prob,
    epiworld_double prob_transmission,
    epiworld_double prob_death,
    epiworld_double prob_noreinfect
    )
{

    ModelSURV(
        *this,
        vname,
        prevalence,
        efficacy_vax,
        latent_period,
        infect_period,
        prob_symptoms,
        prop_vaccinated,
        prop_vax_redux_transm,
        prop_vax_redux_infect,
        surveillance_prob,
        prob_transmission,
        prob_death,
        prob_noreinfect
    );

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/surveillance.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/sirconnected.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SIRCONNECTED_HPP 
#define EPIWORLD_MODELS_SIRCONNECTED_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRCONN : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;

public:

    ModelSIRCONN() {
        
        tracked_agents_infected.reserve(1e4);
        tracked_agents_infected_next.reserve(1e4);

    };

    ModelSIRCONN(
        ModelSIRCONN<TSeq> & model,
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double prob_transmission,
        epiworld_double prob_recovery
    );

    ModelSIRCONN(
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double prob_transmission,
        epiworld_double prob_recovery
    );

    // Tracking who is infected and who is not
    std::vector< epiworld::Agent<TSeq>* > tracked_agents_infected = {};
    std::vector< epiworld::Agent<TSeq>* > tracked_agents_infected_next = {};
    std::vector< epiworld_double >        tracked_agents_weight        = {};
    std::vector< epiworld_double >        tracked_agents_weight_next   = {};

    bool tracked_started = false;
    int tracked_ninfected = 0;
    int tracked_ninfected_next = 0;
    epiworld_double tracked_current_infect_prob = 0.0;

    void run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    Model<TSeq> * clone_ptr();

};

template<typename TSeq>
inline void ModelSIRCONN<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    tracked_agents_infected.clear();
    tracked_agents_infected_next.clear();

    tracked_started = false;
    tracked_ninfected = 0;
    tracked_ninfected_next = 0;
    tracked_current_infect_prob = 0.0;

    Model<TSeq>::run(ndays, seed);

}

template<typename TSeq>
inline Model<TSeq> * ModelSIRCONN<TSeq>::clone_ptr()
{
    
    ModelSIRCONN<TSeq> * ptr = new ModelSIRCONN<TSeq>(
        *dynamic_cast<const ModelSIRCONN<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param prob_transmission Probability of transmission
 * @param prob_recovery Probability of recovery
 */
template<typename TSeq>
inline ModelSIRCONN<TSeq>::ModelSIRCONN(
    ModelSIRCONN<TSeq> & model,
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double prob_transmission,
    epiworld_double prob_recovery
    // epiworld_double prob_reinfection
    )
{


    std::function<void(ModelSIRCONN<TSeq> * m)> tracked_agents_check_init = [](
        ModelSIRCONN<TSeq> * m
        ) -> void
        {

            /* Checking first if it hasn't  */ 
            if (m->tracked_started)
                return;
    
            /* Listing who is infected */ 
            for (auto & p : m->get_agents())
            {
                if (p.get_state() == ModelSIRCONN<TSeq>::INFECTED)
                {
                
                    m->tracked_agents_infected.push_back(&p);
                    m->tracked_ninfected++;
                
                }
            }

            for (auto & p: m->tracked_agents_infected)
            {
                if (p->get_n_viruses() == 0)
                    throw std::logic_error("Cannot be infected and have no viruses.");
            }
            
            m->tracked_started = true;

            // Computing infection probability
            m->tracked_current_infect_prob =  1.0 - std::pow(
                1.0 - (m->par("Contact rate")) * (m->par("Prob. Transmission")) / m->size(),
                m->tracked_ninfected
            );
             

        };

    epiworld::UpdateFun<TSeq> update_susceptible = [
        tracked_agents_check_init
    ](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRCONN<TSeq> * _m = dynamic_cast<ModelSIRCONN<TSeq>*>(m);

            tracked_agents_check_init(_m);

            // No infected individual?
            if (_m->tracked_ninfected == 0)
                return;

            if (m->runif() < _m->tracked_current_infect_prob)
            {

                // Adding the individual to the queue
                _m->tracked_agents_infected_next.push_back(p);
                _m->tracked_ninfected_next++;

                // Now selecting who is transmitting the disease
                epiworld_fast_uint which = static_cast<epiworld_fast_uint>(
                    std::floor(_m->tracked_ninfected * m->runif())
                );


                /* There is a bug in which runif() returns 1.0. It is rare, but
                 * we saw it here. See the Notes section in the C++ manual
                 * https://en.cppreference.com/mwiki/index.php?title=cpp/numeric/random/uniform_real_distribution&oldid=133329
                 * And the reported bug in GCC:
                 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
                 * 
                 */
                if (which == static_cast<epiworld_fast_uint>(_m->tracked_ninfected))
                    --which;

                // Infecting the individual
                p->add_virus(
                    _m->tracked_agents_infected[which]->get_virus(0u),
                    m
                    ); 

                return;

            }

            return;

        };

    epiworld::UpdateFun<TSeq> update_infected = [
        tracked_agents_check_init
    ](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRCONN<TSeq> * _m = dynamic_cast<ModelSIRCONN<TSeq>*>(m);

            tracked_agents_check_init(_m);

            // Is recovering
            if (m->runif() < (m->par("Prob. Recovery")))
            {

                --_m->tracked_ninfected_next;
                epiworld::VirusPtr<int> v = p->get_virus(0u);
                p->rm_virus(0, m);
                return;

            }

            // Will be present next
            _m->tracked_agents_infected_next.push_back(p);

            return;

        };

    epiworld::GlobalFun<TSeq> global_accounting = [](epiworld::Model<TSeq> * m) -> void
        {

            // Getting the right type
            ModelSIRCONN<TSeq> * _m = dynamic_cast<ModelSIRCONN<TSeq>*>(m);

            // On the last day, also reset tracked agents and
            // set the initialized value to false
            if (static_cast<epiworld_fast_uint>(m->today()) == (m->get_ndays() - 1))
            {

                _m->tracked_started = false;
                _m->tracked_agents_infected.clear();
                _m->tracked_agents_infected_next.clear();
                _m->tracked_ninfected = 0;
                _m->tracked_ninfected_next = 0;    
                _m->tracked_current_infect_prob = 0.0;

                return;
            }

            std::swap(_m->tracked_agents_infected, _m->tracked_agents_infected_next);
            _m->tracked_agents_infected_next.clear();

            _m->tracked_ninfected += _m->tracked_ninfected_next;
            _m->tracked_ninfected_next = 0;

            _m->tracked_current_infect_prob = 1.0 - std::pow(
                1.0 - (m->par("Contact rate")) * (m->par("Prob. Transmission")) / m->size(),
                _m->tracked_ninfected
                );

        };

    // Status
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(prob_transmission, "Prob. Transmission");
    model.add_param(prob_recovery, "Prob. Recovery");
    // model.add_param(prob_reinfection, "Prob. Reinfection");
    
    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1, 2, 2);

    model.add_virus(virus, prevalence);

    // Adding updating function
    model.add_global_action(global_accounting, -1);

    model.queuing_off(); // No queuing need

    model.agents_empty_graph(n);

    model.set_name("Susceptible-Infected-Removed (SIR) (connected)");

    return;

}

template<typename TSeq>
inline ModelSIRCONN<TSeq>::ModelSIRCONN(
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double prob_transmission,
    epiworld_double prob_recovery
    )
{

    ModelSIRCONN(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        prob_transmission,
        prob_recovery
    );

    return;

}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/sirconnected.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/seirconnected.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SEIRCONNECTED_HPP
#define EPIWORLD_MODELS_SEIRCONNECTED_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRCONN : public epiworld::Model<TSeq> 
{
public:

    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int RECOVERED   = 3;


    ModelSEIRCONN() {

        tracked_agents_infected.reserve(1e4);
        tracked_agents_infected_next.reserve(1e4);
        
    };

    ModelSEIRCONN(
        ModelSEIRCONN<TSeq> & model,
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double prob_transmission,
        epiworld_double incubation_days,
        epiworld_double prob_recovery
    );
    
    ModelSEIRCONN(
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double prob_transmission,
        epiworld_double incubation_days,
        epiworld_double prob_recovery
    );

    // Tracking who is infected and who is not
    std::vector< epiworld::Agent<>* > tracked_agents_infected = {};
    std::vector< epiworld::Agent<>* > tracked_agents_infected_next = {};

    bool tracked_started = false;
    int tracked_ninfected = 0;
    int tracked_ninfected_next = 0;

    void run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    Model<TSeq> * clone_ptr();

};

template<typename TSeq>
inline void ModelSEIRCONN<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    tracked_agents_infected.clear();
    tracked_agents_infected_next.clear();

    tracked_started = false;
    tracked_ninfected = 0;
    tracked_ninfected_next = 0;

    Model<TSeq>::run(ndays, seed);

}

template<typename TSeq>
inline Model<TSeq> * ModelSEIRCONN<TSeq>::clone_ptr()
{
    
    ModelSEIRCONN<TSeq> * ptr = new ModelSEIRCONN<TSeq>(
        *dynamic_cast<const ModelSEIRCONN<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param prob_transmission Probability of transmission
 * @param prob_recovery Probability of recovery
 */
template<typename TSeq>
inline ModelSEIRCONN<TSeq>::ModelSEIRCONN(
    ModelSEIRCONN<TSeq> & model,
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double prob_transmission,
    epiworld_double incubation_days,
    epiworld_double prob_recovery
    // epiworld_double prob_reinfection
    )
{

    std::function<void(ModelSEIRCONN<TSeq> *)> tracked_agents_check_init = 
    [](ModelSEIRCONN<TSeq> * m) 
        {

            /* Checking first if it hasn't  */ 
            if (m->tracked_started)
                return;

            /* Listing who is infected */ 
            for (auto & p : m->get_agents())
            {
                if (p.get_state() == ModelSEIRCONN<TSeq>::INFECTED)
                {
                
                    m->tracked_agents_infected.push_back(&p);
                    m->tracked_ninfected++;
                
                }
            }

            for (auto & p: m->tracked_agents_infected)
            {
                if (p->get_n_viruses() == 0)
                    throw std::logic_error("Cannot be infected and have no viruses.");
            }
            
            m->tracked_started = true;
                
        };

    epiworld::UpdateFun<TSeq> update_susceptible = 
    [tracked_agents_check_init](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
        {

            // Getting the right type
            ModelSEIRCONN<TSeq> * _m = dynamic_cast<ModelSEIRCONN<TSeq>*>(m);

            tracked_agents_check_init(_m);

            // No infected individual?
            if (_m->tracked_ninfected == 0)
                return;

            // Computing probability of contagion
            // P(infected) = 1 - (1 - beta/Pop * ptransmit) ^ ninfected
            epiworld_double prob_infect = 1.0 - std::pow(
                1.0 - (m->par("Contact rate")) * (m->par("Prob. Transmission")) / m->size(),
                _m->tracked_ninfected
                );

            if (m->runif() < prob_infect)
            {

                // Now selecting who is transmitting the disease
                epiworld_fast_uint which = static_cast<epiworld_fast_uint>(
                    std::floor(_m->tracked_ninfected * m->runif())
                );

                /* There is a bug in which runif() returns 1.0. It is rare, but
                 * we saw it here. See the Notes section in the C++ manual
                 * https://en.cppreference.com/mwiki/index.php?title=cpp/numeric/random/uniform_real_distribution&oldid=133329
                 * And the reported bug in GCC:
                 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
                 * 
                 */
                if (which == static_cast<epiworld_fast_uint>(_m->tracked_ninfected))
                    --which;

                // Infecting the individual
                #ifdef EPI_DEBUG
                if (_m->tracked_agents_infected[which]->get_n_viruses() == 0)
                {

                    printf_epiworld("[epi-debug] date: %i\n", m->today());
                    printf_epiworld("[epi-debug] sim#: %i\n", m->get_n_replicates());

                    throw std::logic_error(
                        "[epi-debug] The agent " + std::to_string(which) + " has no "+
                        "virus to share. The agent's status is: " +
                        std::to_string(_m->tracked_agents_infected[which]->get_state())
                    );
                }
                #endif
                p->add_virus(
                    _m->tracked_agents_infected[which]->get_virus(0u),
                    m,
                    ModelSEIRCONN<TSeq>::EXPOSED
                    ); 

                return;

            }

            return;

        };

    epiworld::UpdateFun<TSeq> update_infected = 
    [tracked_agents_check_init](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
        {

            // Getting the right type
            ModelSEIRCONN<TSeq> * _m = dynamic_cast<ModelSEIRCONN<TSeq>*>(m);

            tracked_agents_check_init(_m);

            auto status = p->get_state();

            if (status == ModelSEIRCONN<TSeq>::EXPOSED)
            {

                // Does the agent become infected?
                if (m->runif() < 1.0/(m->par("Avg. Incubation days")))
                {
                    // Adding the individual to the queue
                    _m->tracked_agents_infected_next.push_back(p);
                    _m->tracked_ninfected_next++;

                    p->change_state(m, ModelSEIRCONN<TSeq>::INFECTED);

                    return;

                }


            } else if (status == ModelSEIRCONN<TSeq>::INFECTED)
            {

                if (m->runif() < (m->par("Prob. Recovery")))
                {

                    _m->tracked_ninfected_next--;
                    p->rm_virus(0, m);
                    return;

                }

                _m->tracked_agents_infected_next.push_back(p);

            } 

            return;

        };

    epiworld::GlobalFun<TSeq> global_accounting = 
    [](epiworld::Model<TSeq>* m) -> void
        {

            // Getting the right type
            ModelSEIRCONN<TSeq> * _m = dynamic_cast<ModelSEIRCONN<TSeq>*>(m);

            // On the last day, also reset tracked agents and
            // set the initialized value to false
            if (static_cast<epiworld_fast_uint>(m->today()) == (m->get_ndays() - 1))
            {

                _m->tracked_started = false;
                _m->tracked_agents_infected.clear();
                _m->tracked_agents_infected_next.clear();
                _m->tracked_ninfected = 0;
                _m->tracked_ninfected_next = 0;    

                return;
            }

            std::swap(_m->tracked_agents_infected, _m->tracked_agents_infected_next);
            _m->tracked_agents_infected_next.clear();

            _m->tracked_ninfected += _m->tracked_ninfected_next;
            _m->tracked_ninfected_next = 0;

        };

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(prob_transmission, "Prob. Transmission");
    model.add_param(prob_recovery, "Prob. Recovery");
    model.add_param(incubation_days, "Avg. Incubation days");
    
    // Status
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Exposed", update_infected);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");


    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(ModelSEIRCONN<TSeq>::EXPOSED, ModelSEIRCONN<TSeq>::RECOVERED, ModelSEIRCONN<TSeq>::RECOVERED);
    model.add_virus(virus, prevalence);

    // Adding updating function
    model.add_global_action(global_accounting, -1);

    model.queuing_off(); // No queuing need

    // Adding the empty population
    model.agents_empty_graph(n);

    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR) (connected)");

    return;

}

template<typename TSeq>
inline ModelSEIRCONN<TSeq>::ModelSEIRCONN(
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double prob_transmission,
    epiworld_double incubation_days,
    epiworld_double prob_recovery
    )
{

    ModelSEIRCONN(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        prob_transmission,
        incubation_days,
        prob_recovery
    );

    return;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/seirconnected.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/seird.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef EPIWORLD_MODELS_SEIRD_HPP
#define EPIWORLD_MODELS_SEIRD_HPP

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed-Deceased (SEIRD) model
*/
template<typename TSeq = int>
class ModelSEIRD : public epiworld::Model<TSeq>
{

protected:

    enum S {
        Susceptible,
        Exposed,
        Infected,
        Hospitalized,
        Recovered,
        Deceased
    };
    
    static void update_exposed(epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m);
    static void update_infected(epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m);
    static void update_hospitalized(epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m);
    static void contact(Model<TSeq> * m);

public:

    ModelSEIRD() {};

    ModelSEIRD(
        ModelSEIRD<TSeq> & model,
        std::string vname
    );

    ModelSEIRD(
        std::string vname,
        epiworld_double prevalence,
        epiworld_double incu_shape,
        epiworld_double incu_rate,
        epiworld_double infe_shape,
        epiworld_double infe_rate,
        epiworld_double p_hosp,
        epiworld_double p_hosp_rec,
        epiworld_double p_hosp_die,
        epiworld_double p_transmission,
        epiworld_double p_transmission_entity,
        size_t n_entities,
        size_t n_interactions
    );

    ModelSEIRD(
        std::string fn,
        std::string vname
        );
    

};

template<typename TSeq>
inline ModelSEIRD<TSeq>::ModelSEIRD(
    std::string fn,
    std::string vname
    )
{

    std::vector<std::string> expected = {
        "Gamma shape (incubation)",
        "Gamma rate (incubation)",
        "Gamma shape (infected)",
        "Gamma rate (infected)",
        "Hospitalization prob.",
        "Prob. hosp. recovers",
        "Prob. hosp. dies",
        "Infectiousness",
        "Infectiousness in entity",
        "Prevalence",
        "N entities",
        "N interactions"
        };

    // Reading the parameters in
    this->read_params(fn);

    for (const auto & p : expected)
    {
        if (this->parameters.find(p) == this->parameters.end())
        {
            throw std::logic_error(
                std::string(
                    "The specified file doesn't have all the required parameters for the model"
                )
            );
        }
    }


    // Setting up the model
    ModelSEIRD(*this, vname);

    return;

}

template<typename TSeq>
inline ModelSEIRD<TSeq>::ModelSEIRD(
    ModelSEIRD<TSeq> & model,
    std::string vname
    )
{

    model.add_state(
        "Susceptible", 
        sampler::make_update_susceptible<TSeq>({S::Exposed, S::Hospitalized})
        );

    model.add_state("Exposed", this->update_exposed);
    model.add_state("Infected", this->update_infected);
    model.add_state("Hospitalized", this->update_hospitalized);
    model.add_state("Recovered");
    model.add_state("Deceased");

    // Creating the virus
    Virus<TSeq> virus(vname);
    virus.set_prob_infecting(&model("Infectiousness"));
    virus.set_state(S::Exposed, S::Recovered);
    virus.set_queue(QueueValues::OnlySelf, -99LL);
    virus.get_data() = {0.0, 0.0F};

    model.add_virus_n(virus, model("Prevalence"));
     
    // Adding randomly distributed entities, each one with capacity for entity_capacity
    for (size_t r = 0u; r < model("N entities"); ++r)
    {
        Entity<TSeq> e(std::string("Location ") + std::to_string(r));
        model.add_entity(e);
    }

    // This will act through the global
    model.add_global_action(contact, -99);

    model.set_name("Susceptible-Exposed-Infected-Recovered-Deceased (SEIRD)");

    return;
   
}

template<typename TSeq>
inline ModelSEIRD<TSeq>::ModelSEIRD(
    std::string vname,
    epiworld_double prevalence,
    epiworld_double incu_shape,
    epiworld_double incu_rate,
    epiworld_double infe_shape,
    epiworld_double infe_rate,
    epiworld_double p_hosp,
    epiworld_double p_hosp_rec,
    epiworld_double p_hosp_die,
    epiworld_double p_transmission,
    epiworld_double p_transmission_entity,
    size_t n_entities,
    size_t n_interactions
) {

    // Adding parameters
    this->add_param(incu_shape, "Gamma shape (incubation)");
    this->add_param(incu_rate, "Gamma rate (incubation) ");
    this->add_param(infe_shape, "Gamma shape (infected)");
    this->add_param(infe_rate, "Gamma rate (infected)");
    this->add_param(p_hosp, "Hospitalization prob.");
    this->add_param(p_hosp_rec, "Prob. hosp. recovers");
    this->add_param(p_hosp_rec, "Prob. hosp. dies");
    this->add_param(p_transmission, "Infectiousness");
    this->add_param(p_transmission_entity, "Infectiousness in entity");
    this->add_param(prevalence, "Prevalence");
    this->add_param(static_cast<epiworld_double>(n_entities), "N entities");
    this->add_param(static_cast<epiworld_double>(n_interactions), "N interactions");

    ModelSEIRD<TSeq>(
        *this,
        vname
        );

    return;

}




// Update dynamics for exposed individuals
template<typename TSeq>
inline void ModelSEIRD<TSeq>::update_exposed(
    epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
)
{
    // Checking if the data was assigned
    auto v = p->get_virus(0u);
    if (m->today() == (v->get_date() + 1))
    {

        // Checking the same
        if (v->get_data().size() != 3u)
            v->get_data().resize(3u);

        // Days of incubation
        v->get_data()[0u] = v->get_date() + m->rgamma(
            m->par("Gamma shape (incubation)"),
            m->par("Gamma rate (incubation)")
            );

        // Duration as Infected
        v->get_data()[1u] = v->get_date() + m->rgamma(
            m->par("Gamma shape (infected)"),
            m->par("Gamma rate (infected)")
            );

        // Prob of becoming hospitalized
        v->get_data()[2u] = (
            m->runif() < m->par("Hospitalization prob.")
            ) ? 100.0 : -100.0;


    } else if (m->today() >= v->get_data()[0u])
    {
        p->change_state(m, S::Infected, epiworld::QueueValues::Everyone);
        return;
    }

    return;

}

// Update dynamics for infected individuals
template<typename TSeq>
inline void ModelSEIRD<TSeq>::update_infected(
    epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
)
{

    // Computing probability of recovery
    auto v = p->get_virus(0u);

    if (m->today() < v->get_data()[1u])
        return;

    if (v->get_data()[2u] < 0)
    {
        
        p->rm_virus(v, m, S::Recovered, -epiworld::QueueValues::Everyone);
        return;

    }
    else
    {

        // Individual goes hospitalized
        p->change_state(m, S::Hospitalized, -epiworld::QueueValues::Everyone);
        return;

    }

}

// Update dynamics for hospitalized individuals
template<typename TSeq>
inline void ModelSEIRD<TSeq>::update_hospitalized(
    epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
)
{

    // Computing the recovery probability
    auto v = p->get_virus(0u);
    auto probs = {
        m->par("Prob. hosp. recovers"),
        m->par("Prob. hosp. dies")
        };

    int which = epiworld::roulette(probs, m);

    // Nothing happens
     if (which < 0)
        return;

    if (which == 0) // Then it recovered
    {
        p->rm_virus(v, m, S::Recovered, epiworld::QueueValues::NoOne);
        return;
    }

    // Individual dies
    p->rm_virus(v, m, S::Deceased, epiworld::QueueValues::NoOne);

    return;

}

/**
 * @brief Transmission by contact outside home
 */
template<typename TSeq>
inline void ModelSEIRD<TSeq>::contact(Model<TSeq> * m)
{
    for (auto & a : (m->get_agents()))
    {
        // Will it get it from the entities?
        if (a.get_state() == S::Susceptible)
        {

            AgentsSample<int> neighbors(m, a, m->par("N interactions"), true);

            int n_viruses = 0;
            for (auto n : neighbors) {
                if (n->get_state() == S::Infected)
                    m->array_virus_tmp[n_viruses++] = &(*n->get_virus(0u));
            }

            // Nothing to see here
            if (n_viruses == 0)
                continue;

            // Is the individual getting the infection?
            double p_infection = 1.0 - std::pow(
	        1.0 - m->par("Infectiousness in entity"), n_viruses
		);

            if (m->runif() >= p_infection)
                continue;

            // Who infects the individual?
            int which = std::floor(m->runif() * n_viruses);

            a.add_virus(
                *(m->array_virus_tmp[which]), // Viruse.
                m, 
                S::Exposed,                   // New state.
                QueueValues::OnlySelf         // Change on the queue.
                ); 

        }


    }
}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/seird.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/epiworld//models/sirlogit.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "../epiworld.hpp"

#ifndef EPIWORLD_MODELS_SIRLOGIT_HPP 
#define EPIWORLD_MODELS_SIRLOGIT_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRLogit : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;

public:

    ModelSIRLogit() {};

    /**
      * @param vname Name of the virus.
      * @param coefs_infect Double ptr. Infection coefficients.
      * @param coefs_recover Double ptr. Recovery coefficients.
      * @param ncoef_infect Unsigned int. Number of infection coefficients.
      * @param ncoef_recover Unsigned int. Number of recovery coefficients.
      * @param coef_infect_cols Vector<unsigned int>. Ids of infection vars.
      * @param coef_recover_cols Vector<unsigned int>. Ids of recover vars.
    */
    ModelSIRLogit(
        ModelSIRLogit<TSeq> & model,
        std::string vname,
        double * data,
        size_t ncols,
        std::vector< double > coefs_infect,
        std::vector< double > coefs_recover,
        std::vector< size_t > coef_infect_cols,
        std::vector< size_t > coef_recover_cols,
        epiworld_double prevalence
    );

    ModelSIRLogit(
        std::string vname,
        double * data,
        size_t ncols,
        std::vector< double > coefs_infect,
        std::vector< double > coefs_recover,
        std::vector< size_t > coef_infect_cols,
        std::vector< size_t > coef_recover_cols,
        epiworld_double prevalence
    );

    void run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    Model<TSeq> * clone_ptr();

    bool tracked_started = false;
    
    std::vector< double > coefs_infect;
    std::vector< double > coefs_recover;
    std::vector< size_t > coef_infect_cols;
    std::vector< size_t > coef_recover_cols;

};



template<typename TSeq>
inline void ModelSIRLogit<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    Model<TSeq>::run(ndays, seed);

}

template<typename TSeq>
inline Model<TSeq> * ModelSIRLogit<TSeq>::clone_ptr()
{
    
    ModelSIRLogit<TSeq> * ptr = new ModelSIRLogit<TSeq>(
        *dynamic_cast<const ModelSIRLogit<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param prob_transmission Probability of transmission
 * @param prob_recovery Probability of recovery
 */
template<typename TSeq>
inline ModelSIRLogit<TSeq>::ModelSIRLogit(
    ModelSIRLogit<TSeq> & model,
    std::string vname,
    double * data,
    size_t ncols,
    std::vector< double > coefs_infect,
    std::vector< double > coefs_recover,
    std::vector< size_t > coef_infect_cols,
    std::vector< size_t > coef_recover_cols,
    epiworld_double prevalence
    )
{

    if (coef_infect_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_infect_cols.");

    if (coef_recover_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_recover_cols.");

    // Saving the variables
    model.set_agents_data(
        data, ncols
    );

    model.coefs_infect = coefs_infect;
    model.coefs_recover = coefs_recover;
    model.coef_infect_cols = coef_infect_cols;
    model.coef_recover_cols = coef_recover_cols;

    std::function<void(ModelSIRLogit<TSeq> * m)> check_init = [](
        ModelSIRLogit<TSeq> * m
        ) -> void
        {

            /* Checking first if it hasn't  */ 
            if (m->tracked_started)
                return;

            /* Checking specified columns in the model */
            for (const auto & c : m->coef_infect_cols)
            {
                if (c >= m->agents_data_ncols)
                    throw std::range_error("Columns specified in coef_infect_cols out of range.");
            }

            for (const auto & c : m->coef_recover_cols)
            {
                if (c >= m->agents_data_ncols)
                    throw std::range_error("Columns specified in coef_recover_cols out of range.");
            }
    
            /* Checking attributes */ 
            if (m->coefs_infect.size() != (m->coef_infect_cols.size() + 1u))
                throw std::logic_error(
                    "The number of coefficients (infection) doesn't match the number of features. It must be as many features of the agents plus 1 (exposure.)"
                    );

            if (m->coefs_recover.size() != m->coef_recover_cols.size())
                throw std::logic_error(
                    "The number of coefficients (recovery) doesn't match the number of features. It must be as many features of the agents."
                    );
            
            m->tracked_started = true;            

        };

    epiworld::UpdateFun<TSeq> update_susceptible = [
        check_init
    ](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRLogit<TSeq> * _m = dynamic_cast<ModelSIRLogit<TSeq>*>(m);

            check_init(_m);

            // Exposure coefficient
            const double coef_exposure = _m->coefs_infect[0u];

            // This computes the prob of getting any neighbor variant
            size_t nvariants_tmp = 0u;

            size_t id      = p->get_id();
            size_t nagents = m->size();

            double baseline = 0.0;
            for (size_t k = 0u; k < _m->coef_infect_cols.size(); ++k)
            {
                baseline += (*(
                        m->get_agents_data() +
                        /* data is stored column-major */
                        (_m->coef_infect_cols[k] * nagents + id)
                    )) * _m->coefs_infect[k + 1u];
            }

            for (auto & neighbor: p->get_neighbors()) 
            {
                
                for (const VirusPtr<TSeq> & v : neighbor->get_viruses()) 
                { 

                    #ifdef EPI_DEBUG
                    if (nvariants_tmp >= m->array_virus_tmp.size())
                        throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                    #endif
                        
                    /* And it is a function of susceptibility_reduction as well */ 
                    m->array_double_tmp[nvariants_tmp] =
                        baseline +
                        (1.0 - p->get_susceptibility_reduction(v, m)) * 
                        v->get_prob_infecting(m) * 
                        (1.0 - neighbor->get_transmission_reduction(v, m))  *
                        coef_exposure
                        ; 

                    // Applying the plogis function
                    m->array_double_tmp[nvariants_tmp] = 1.0/
                        (1.0 + std::exp(-m->array_double_tmp[nvariants_tmp]));
                
                    m->array_virus_tmp[nvariants_tmp++] = &(*v);

                }

            }

            // No virus to compute
            if (nvariants_tmp == 0u)
                return;

            // Running the roulette
            int which = roulette(nvariants_tmp, m);

            if (which < 0)
                return;

            p->add_virus(*m->array_virus_tmp[which], m);

            return;

        };

    epiworld::UpdateFun<TSeq> update_infected = [
        check_init
    ](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRLogit<TSeq> * _m = dynamic_cast<ModelSIRLogit<TSeq>*>(m);

            check_init(_m);

            // Computing recovery probability once
            double prob    = 0.0;
            size_t id      = p->get_id();
            size_t nagents = m->size();
            #pragma omp simd reduction(+:prob)
            for (size_t i = 0u; i < _m->coefs_recover.size(); ++i)
            {
                prob +=
                    (*(
                        m->get_agents_data() +
                        /* data is stored column-major */
                        (_m->coef_recover_cols[i] * nagents + id)
                    )) * _m->coefs_recover[i];

            }

            // Computing logis
            prob = 1.0/(1.0 + std::exp(-prob));

            if (prob > m->runif())
                p->rm_virus(0, m);
            
            return;

        };

    // Status
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Setting up parameters
    // model.add_param(contact_rate, "Contact rate");
    // model.add_param(prob_transmission, "Prob. Transmission");
    // model.add_param(prob_recovery, "Prob. Recovery");
    // model.add_param(prob_reinfection, "Prob. Reinfection");
    
    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(
        ModelSIRLogit<TSeq>::INFECTED,
        ModelSIRLogit<TSeq>::RECOVERED,
        ModelSIRLogit<TSeq>::RECOVERED
        );

    // virus.set_prob

    model.add_virus(virus, prevalence);

    model.set_name("Susceptible-Infected-Removed (SIR) (logit)");

    return;

}

template<typename TSeq>
inline ModelSIRLogit<TSeq>::ModelSIRLogit(
    std::string vname,
    double * data,
    size_t ncols,
    std::vector< double > coefs_infect,
    std::vector< double > coefs_recover,
    std::vector< size_t > coef_infect_cols,
    std::vector< size_t > coef_recover_cols,
    epiworld_double prevalence
    )
{

    ModelSIRLogit(
        *this,
        vname,
        data,
        ncols,
        coefs_infect,
        coefs_recover,
        coef_infect_cols,
        coef_recover_cols,
        prevalence
    );

    return;

}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld//models/sirlogit.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/epiworld/models/models.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



}

#endif 
