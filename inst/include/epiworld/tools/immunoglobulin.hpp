#ifndef EPIWORLD_TOOLS_IMMUNOGLOBULIN_HPP
#define EPIWORLD_TOOLS_IMMUNOGLOBULIN_HPP

#include <memory>
#include <vector>
#include <cassert>
#include <string_view>
#include "../config.hpp"
#include "../tool-bones.hpp"

/**
 * @brief Template for a common Immunoglobulin tool.
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
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ToolImmunoglobulin: public Tool<TSeq> {
private:

    int _immune = -1;
    std::string _par_efficacy;
    epiworld_double _set_immunity(Model<TSeq> & model);

    // This half-life
    std::string _par_half_life_mean;
    std::string _par_half_life_sd;
    int _removal_time = -1;
    bool _remove_tool();

    // This method is used to throw an error when trying to set a method
    // that cannot be set for this tool.
    void _error_method(std::string_view method_name) const;

public:
    ToolImmunoglobulin(
        std::string name,
        std::string par_efficacy,
        std::string par_half_life_mean,
        std::string par_half_life_sd
    ) : Tool<TSeq>(name) {
        this->_par_efficacy = par_efficacy;
        this->_par_half_life_mean = par_half_life_mean;
        this->_par_half_life_sd = par_half_life_sd;
    };

    virtual epiworld_double get_susceptibility_reduction(
        VirusPtr<TSeq> & v,
        Model<TSeq> * model
    ) override;
    
    virtual void set_susceptibility_reduction_fun(ToolFun<TSeq> fun) override;
    virtual void set_susceptibility_reduction(std::string param) override;
    virtual void set_susceptibility_reduction(epiworld_double prob) override;

    std::unique_ptr<Tool<TSeq>> clone_ptr() const override;

};

template<typename TSeq>
inline void ToolImmunoglobulin<TSeq>::_error_method(std::string_view method_name) const {
    throw std::logic_error(
        std::string(
            "The method " + std::string(method_name) + " cannot be set for the "
        ) +
        std::string("ToolImmunoglobulin (Tool).")
    );
}

template<typename TSeq>
inline epiworld_double ToolImmunoglobulin<TSeq>::get_susceptibility_reduction(
    VirusPtr<TSeq> &,
    Model<TSeq> * model
)
{

    // Updating a single agent (if needed)
    if (_immune == -1)
    {
        _set_immunity(*model);
    }

    // Deciding if the tool should be removed
    #ifdef EPI_DEBUG
    EPI_DEBUG_PRINTF(
        "Agent %d; sim_id %li; today %d; removal_time %d\n",
        this->agent->get_id(), model->get_sim_id(),
        model->today(), _removal_time
    );
    #endif
    if (model->today() >= _removal_time)
    {
        this->get_agent()->rm_tool(*model, this->pos_in_agent);
        return 0.0;
    }

    return  (_immune == 1) ? 1.0 : 0.0;

}

template<typename TSeq>
inline void ToolImmunoglobulin<TSeq>::set_susceptibility_reduction_fun(
    ToolFun<TSeq>
)
{
    _error_method("set_susceptibility_reduction_fun");
}

template<typename TSeq>
inline void ToolImmunoglobulin<TSeq>::set_susceptibility_reduction(std::string)
{
    _error_method("set_susceptibility_reduction(string)");
}

template<typename TSeq>
inline void ToolImmunoglobulin<TSeq>::set_susceptibility_reduction(epiworld_double)
{

    _error_method("set_susceptibility_reduction(epiworld_double)");

}

template<typename TSeq>
inline epiworld_double ToolImmunoglobulin<TSeq>::_set_immunity(Model<TSeq> & model)
{

    // Setting immunity
    _immune = (model.runif() < model.par(this->_par_efficacy)) ? 1 : 0;

    // Setting the removal time (regardless)
    auto half_life = model.par(this->_par_half_life_mean);
    if (half_life > 0.0)
    {

        // Drawing the half-life for this agent
        auto hl_sd = model.par(this->_par_half_life_sd);
        auto hl_mean = model.par(this->_par_half_life_mean);

        if (hl_sd > 0.0)
            half_life = model.rnorm(hl_mean, hl_sd);
        else
            half_life = hl_mean;

        _removal_time = model.today() + static_cast<int>(std::round(half_life));
    }

    return _immune;
}

template<typename TSeq>
inline std::unique_ptr<Tool<TSeq>> ToolImmunoglobulin<TSeq>::clone_ptr() const
{
    auto ans =  std::make_unique<ToolImmunoglobulin<TSeq>>(*this);
    ans->_immune = -1;
    return ans;
}

#endif