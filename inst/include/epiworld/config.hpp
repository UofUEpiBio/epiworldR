#ifndef EPIWORLD_CONFIG_HPP
#define EPIWORLD_CONFIG_HPP

#ifndef printf_epiworld
    #define printf_epiworld fflush(stdout);printf
#endif

// In case the user has a way to stop the program
// This is called during `run_multiple()` and it is
// passed the simulation number.
#ifndef EPI_CHECK_USER_INTERRUPT
    #define EPI_CHECK_USER_INTERRUPT(a)
#endif

#ifndef EPIWORLD_MAXNEIGHBORS
    #define EPIWORLD_MAXNEIGHBORS 1048576
#endif

#if defined(_OPENMP) || defined(__OPENMP)
    #include <omp.h>
// #else
//     #define omp_get_thread_num() 0
//     #define omp_set_num_threads() 1
#endif

#ifndef epiworld_double
    #define epiworld_double float
#endif

#ifndef epiworld_fast_int
    #define epiworld_fast_int int
#endif

#ifndef epiworld_fast_uint
    #define epiworld_fast_uint unsigned long long int
#endif

#define EPI_DEFAULT_TSEQ int

#ifndef EPI_MAX_TRACKING
    #define EPI_MAX_TRACKING 100
#endif

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Model;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Agent;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class PersonTools;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Virus;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Viruses;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Viruses_const;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Tool;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Tools;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Tools_const;

template<typename TSeq = EPI_DEFAULT_TSEQ>
class Entity;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using VirusPtr = std::shared_ptr< Virus< TSeq > >;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using ToolPtr = std::shared_ptr< Tool< TSeq > >;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using ToolFun = std::function<epiworld_double(Tool<TSeq>&,Agent<TSeq>*,VirusPtr<TSeq>,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using MixerFun = std::function<epiworld_double(Agent<TSeq>*,VirusPtr<TSeq>,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using MutFun = std::function<bool(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using PostRecoveryFun = std::function<void(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using VirusFun = std::function<epiworld_double(Agent<TSeq>*,Virus<TSeq>&,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using UpdateFun = std::function<void(Agent<TSeq>*,Model<TSeq>*)>;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using GlobalFun = std::function<void(Model<TSeq>*)>;

template<typename TSeq>
struct Event;

template<typename TSeq = EPI_DEFAULT_TSEQ>
using EventFun = std::function<void(Event<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute viruses at initialization
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
using VirusToAgentFun = std::function<void(Virus<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute tools at initialization
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
using ToolToAgentFun = std::function<void(Tool<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Decides how to distribute entities at initialization
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
using EntityToAgentFun = std::function<void(Entity<TSeq>&,Model<TSeq>*)>;

/**
 * @brief Event data for update an agent
 * 
 * @tparam TSeq 
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
struct Event {
    Agent<TSeq> * agent;
    VirusPtr<TSeq> virus;
    ToolPtr<TSeq> tool;
    Entity<TSeq> * entity;
    epiworld_fast_int new_state;
    epiworld_fast_int queue;
    EventFun<TSeq> call;
    int idx_agent;
    int idx_object;
public:
/**
     * @brief Construct a new Event object
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
    Event(
        Agent<TSeq> * agent_,
        VirusPtr<TSeq> & virus_,
        ToolPtr<TSeq> & tool_,
        Entity<TSeq> * entity_,
        epiworld_fast_int new_state_,
        epiworld_fast_int queue_,
        EventFun<TSeq> & call_,
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

#ifndef EPI_DEFAULT_INCUBATION_DAYS
    #define EPI_DEFAULT_INCUBATION_DAYS         7.0
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

    #define epiexception(a) std::logic_error
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
    #define epiexception(a) a
#endif

#if defined(EPI_DEBUG_NO_THREAD_ID) || (!defined(__OPENMP) && !defined(_OPENMP))
    #define EPI_GET_THREAD_ID() 0
#else
    #define EPI_GET_THREAD_ID() omp_get_thread_num()
#endif

#endif