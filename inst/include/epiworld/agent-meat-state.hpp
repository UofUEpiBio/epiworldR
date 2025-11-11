#ifndef EPIWORLD_PERSON_MEAT_STATE_HPP
#define EPIWORLD_PERSON_MEAT_STATE_HPP

#include "model-bones.hpp"
#include "agent-meat-virus-sampling.hpp"
#include "config.hpp"

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

template<typename TSeq = EPI_DEFAULT_TSEQ>
inline void default_update_susceptible(
    Agent<TSeq> * p,
    Model<TSeq> * m
    )
{

    Virus<TSeq> * virus = sampler::sample_virus_single<TSeq>(p, m);
    
    if (virus == nullptr)
        return;

    p->set_virus(*virus, m); 

    return;

}

template<typename TSeq = EPI_DEFAULT_TSEQ>
inline void default_update_exposed(Agent<TSeq> * p, Model<TSeq> * m) {

    if (p->get_virus() == nullptr)
        throw std::logic_error(
            std::string("Using the -default_update_exposed- on agents WITHOUT viruses makes no sense! ") +
            std::string("Agent id ") + std::to_string(p->get_id()) + std::string(" has no virus registered.")
            );

    // Die
    auto & virus = p->get_virus();
    m->array_double_tmp[0u] = 
        virus->get_prob_death(m) * (1.0 - p->get_death_reduction(virus, m)); 

    // Recover
    m->array_double_tmp[1u] = 
        1.0 - (1.0 - virus->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(virus, m)); 
    

    // Running the roulette
    int which = roulette(2u, m);

    if (which < 0)
        return;

    // Which roulette happen?
    if (which == 0u) // If odd
    {

        p->rm_agent_by_virus(m);
        
    } else {

        p->rm_virus(m);

    }

    return ;

};

#endif