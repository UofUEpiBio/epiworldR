#ifndef QUARANTINETRIGGER_HPP
#define QUARANTINETRIGGER_HPP

#include <vector>
#include "../config.hpp"
#include "quarantinetrigger-bones.hpp"

template<typename TSeq = EPI_DEFAULT_TSEQ>
class QuarantineTrigger {
private:
    int _model_sim_id = -1;
    int _day = -1;
    std::vector< size_t > _agents_triggering_quarantine;
    std::vector< int > _date_infectious;
    
    void _setup(const Model<TSeq> & model);
public:
    QuarantineTrigger() = default;
    void add_triggering_agent(
        const Model<TSeq> & model,
        const Agent<TSeq> & agent,
        int date_infectious
    );

    std::vector< size_t > & get_triggering_agents();
    std::vector< int > & get_date_infectious();
};

template<typename TSeq>
inline void QuarantineTrigger<TSeq>::_setup(const Model<TSeq> & model) {

    if (
        (static_cast<int>(model.get_sim_id()) != _model_sim_id) ||
        (static_cast<int>(model.today()) != _day)
    ) {
        _model_sim_id = static_cast<int>(model.get_sim_id());
        _day = static_cast<int>(model.today());
        _agents_triggering_quarantine.clear();
        _date_infectious.clear();
    }

}

template<typename TSeq>
inline void QuarantineTrigger<TSeq>::add_triggering_agent(
    const Model<TSeq> & model,
    const Agent<TSeq> & agent,
    int date_infectious
) {
    
    // Ensuring we are not in a different run
    _setup(model);
    this->_agents_triggering_quarantine.push_back(agent.get_id());
    this->_date_infectious.push_back(date_infectious);
}

template<typename TSeq>
inline std::vector< size_t > & QuarantineTrigger<TSeq>::get_triggering_agents() {
    return this->_agents_triggering_quarantine;
}

template<typename TSeq>
inline std::vector< int > & QuarantineTrigger<TSeq>::get_date_infectious() {
    return this->_date_infectious;
}

#endif