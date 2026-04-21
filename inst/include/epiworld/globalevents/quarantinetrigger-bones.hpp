#ifndef QUARANTINETRIGGER_BONES_HPP
#define QUARANTINETRIGGER_BONES_HPP

#include <vector>
#include "../config.hpp"

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

#endif