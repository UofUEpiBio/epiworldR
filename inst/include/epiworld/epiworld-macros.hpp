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


#define EPI_NEW_ENTITYTOAGENTFUN(funname,tseq) inline void \
    (funname)(epiworld::Entity<tseq> & e, epiworld::Model<tseq> * m)

#define EPI_NEW_ENTITYTOAGENTFUN_LAMBDA(funname,tseq) \
    epiworld::EntityToAgentFun<tseq> funname = \
    [](epiworld::Entity<tseq> & e, epiworld::Model<tseq> * m) -> void

// Use this to make it more efficient for storage if the type is small
// and the pointer is not needed
#define EPI_TYPENAME_TRAITS(tseq, bound) typename std::conditional< \
    sizeof( tseq ) <= sizeof( bound ), \
    tseq, \
    std::shared_ptr< tseq > \
    >::type

#define EPI_IF_TSEQ_LESS_EQ_INT(tseq) \
    if constexpr (sizeof( tseq ) <= sizeof( int ))

#define EPI_CHECK_COALESCE(proposed_, virus_tool_, alt_) \
    if (static_cast<int>(proposed_) == -99) {\
        if (static_cast<int>(virus_tool_) == -99) \
            (proposed_) = (alt_);\
        else (proposed_) = (virus_tool_);}

#endif
