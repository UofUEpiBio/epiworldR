#ifndef EPIWORLD_TOOLS_VACCINE_HPP
#define EPIWORLD_TOOLS_VACCINE_HPP

#include <memory>
#include <vector>
#include <cassert>
#include "../config.hpp"
#include "../tool-bones.hpp"

/**
 * @brief Template for a common vaccine
 * 
 * The main difference with a regular tool is that the reduction
 * in susceptibility is at the agent level and fixed for the
 * entire simulation. In other words, if the efficacy is 65%,
 * then 65% of the agents that receive the vaccine will be fully
 * protected (100% reduction in susceptibility) and 35% will not be
 * protected at all (0% reduction in susceptibility).
 * 
 * @ingroup tools
 * 
 */
class ToolVaccine: public Tool {
private:

    int _immune = -1;
    epiworld_double _efficacy = 0.0;
    epiworld_double _set_immunity(Model & model);

public:
    ToolVaccine(std::string name = "Vaccine") : Tool(name) {};

    virtual epiworld_double get_susceptibility_reduction(
        VirusPtr & v,
        Model * model
    ) override;
    
    virtual void set_susceptibility_reduction_fun(ToolFun fun) override;
    virtual void set_susceptibility_reduction(std::string param) override;
    virtual void set_susceptibility_reduction(epiworld_double prob) override;

    std::unique_ptr<Tool> clone_ptr() const override;

};

inline epiworld_double ToolVaccine::get_susceptibility_reduction(
    VirusPtr &,
    Model * model
)
{

    // Updating a single agent (if needed)
    if (_immune == -1)
    {
        _set_immunity(*model);
    }

    return  (_immune == 1) ? 1.0 : 0.0;

}

inline void ToolVaccine::set_susceptibility_reduction_fun(
    ToolFun
)
{
    throw std::logic_error(
        std::string(
            "The susceptibility reduction function cannot be set for the "
         ) +
        std::string("ToolVaccine (Tool).")
    );
}

inline void ToolVaccine::set_susceptibility_reduction(std::string)
{
    throw std::logic_error(
        std::string(
            "The susceptibility reduction probability cannot be set for the "
        ) +
        std::string(
            "ToolVaccine (Tool) using a pointer. Use the version that takes a "
        ) +
        std::string("double instead.")
    );
}

inline void ToolVaccine::set_susceptibility_reduction(epiworld_double prob)
{

    // assertm(prob >= 0 && prob <= 1, "The efficacy must be between 0 and 1.");
    _efficacy = prob;

}

inline epiworld_double ToolVaccine::_set_immunity(Model & model)
{
    _immune = (model.runif() < _efficacy) ? 1 : 0;
    return _immune;
}

inline std::unique_ptr<Tool> ToolVaccine::clone_ptr() const
{
    auto ans =  std::make_unique<ToolVaccine>(*this);
    ans->_immune = -1;
    return ans;
}

#endif