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
#include "agent-meat-virus-sampling.hpp"


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