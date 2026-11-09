#ifndef EPIWORLD_MACROS_HPP
#define EPIWORLD_MACROS_HPP

/**
 * @brief Helper macro to define a new tool
 * 
 */
#define EPI_NEW_TOOL(fname) inline epiworld_double \
(fname)(\
    epiworld::Tool & t, \
    epiworld::Agent * p, \
    std::shared_ptr<epiworld::Virus> v, \
    epiworld::Model * m\
    )

/**
 * @brief Create a Tool within a function
 * 
 */
#define EPI_NEW_TOOL_LAMBDA(funname) \
    epiworld::ToolFun funname = \
    [](epiworld::Tool & t, \
    epiworld::Agent * p, \
    std::shared_ptr<epiworld::Virus> v, \
    epiworld::Model * m) -> epiworld_double

/**
 * @brief Helper macro for accessing model parameters
 * 
 */
#define EPI_PARAMS(i) m->operator()(i)

/**
 * @brief Value of the parameter named `pname` (a string literal) in `model`
 * (a pointer to a model).
 *
 * @details Same value as `model->par(pname)`, but the name is looked up once
 * per model layout and cached in a function-local `static ParamRef`, so
 * repeated calls read the value directly. Use it in update functions and
 * other code that runs for many agents every day.
 */
#define EPI_PAR(model, pname) \
    ([]() -> const ::epiworld::ParamRef & { \
        static const ::epiworld::ParamRef epi_param_ref_(pname); \
        return epi_param_ref_; \
    }()(*(model)))

/**
 * @brief Helper macro for defining Mutation Functions
 * 
 */
#define EPI_NEW_MUTFUN(funname) inline bool \
    (funname)(\
    epiworld::Agent * p, \
    epiworld::Virus & v, \
    epiworld::Model * m )

#define EPI_NEW_MUTFUN_LAMBDA(funname) \
    epiworld::MutFun funname = \
    [](epiworld::Agent * p, \
    epiworld::Virus & v, \
    epiworld::Model * m) -> bool

#define EPI_RUNIF() m->runif()

#define EPI_TOKENPASTE(a,b) a ## b
#define MPAR(num) *(m->EPI_TOKENPASTE(p,num))


#define EPI_NEW_GLOBALFUN(funname) inline void \
    (funname)(epiworld::Model* m)

#define EPI_NEW_GLOBALFUN_LAMBDA(funname) \
    epiworld::GlobalFun funname = \
    [](epiworld::Model* m) -> void


#define EPI_CHECK_COALESCE(proposed_, virus_tool_, alt_) \
    if (static_cast<int>(proposed_) == -99) {\
        if (static_cast<int>(virus_tool_) == -99) \
            (proposed_) = (alt_);\
        else (proposed_) = (virus_tool_);}

#endif
