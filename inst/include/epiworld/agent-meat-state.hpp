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

    p->set_virus(*m, *virus);

    return;

}

/**
 * @brief Factory function to create a state update function based on
 * transition rates.
 *
 * @details
 * This function creates an `UpdateFun` that transitions agents between
 * states using named model parameters as per-step transition probabilities.
 * At each time step the generated function reads the current values of the
 * parameters from the model and uses the roulette sampling algorithm to
 * decide whether a transition occurs (and if so, to which target state).
 *
 * The parameter names must correspond to parameters that have already been
 * added to the model (via `Model::add_param()`). The target states are
 * specified as integer state codes (the indices returned by
 * `Model::add_state()`).
 *
 * @tparam TSeq Type of the sequence (default `EPI_DEFAULT_TSEQ`).
 * @param param_names  Names of the model parameters that hold the
 *                     transition rates. The i-th name gives the rate for
 *                     transitioning to the i-th target state.
 * @param target_states Integer codes of the destination states, in the
 *                      same order as `param_names`.
 * @return An `UpdateFun<TSeq>` suitable for use with
 *         `Model::add_state()`.
 *
 * @throws std::logic_error if `param_names` and `target_states` differ in
 *         size, or if either is empty.
 *
 * ### Example: SEIRD model
 *
 * @code{.cpp}
 * // E -> I transition
 * auto update_exposed = new_state_update_transition<int>(
 *     {"E->I transition rate"},
 *     {2}  // state code for Infected
 * );
 *
 * // I -> R or I -> D transitions
 * auto update_infected = new_state_update_transition<int>(
 *     {"I->R transition rate", "I->D transition rate"},
 *     {3, 4}  // state codes for Recovered and Deceased
 * );
 *
 * model.add_state("Susceptible", default_update_susceptible<int>);
 * model.add_state("Exposed",     update_exposed);
 * model.add_state("Infected",    update_infected);
 * model.add_state("Recovered");
 * model.add_state("Deceased");
 * @endcode
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline UpdateFun<TSeq> new_state_update_transition(
    std::vector< std::string > param_names,
    std::vector< epiworld_fast_uint > target_states
)
{

    if (param_names.size() != target_states.size())
        throw std::logic_error(
            "The number of parameter names (" +
            std::to_string(param_names.size()) +
            ") and target states (" +
            std::to_string(target_states.size()) +
            ") must match."
        );

    if (param_names.empty())
        throw std::logic_error(
            "At least one transition must be specified."
        );

    return [param_names, target_states](
        Agent<TSeq> * p,
        Model<TSeq> * m
    ) -> void {

        size_t n = param_names.size();
        int which;

        if (n <= 1024u)
        {
            for (size_t i = 0u; i < n; ++i)
                m->array_double_tmp[i] = m->par(param_names[i]);

            // Roulette sampling: returns -1 if no transition occurs,
            // otherwise the index of the transition that fires.
            which = roulette(static_cast<epiworld_fast_uint>(n), m);
        }
        else
        {
            std::vector< epiworld_double > probs(n);
            for (size_t i = 0u; i < n; ++i)
                probs[i] = m->par(param_names[i]);

            // Fallback for transition tables larger than the temporary buffer.
            which = roulette(probs, m);
        }
        if (which < 0)
            return;

        p->change_state(*m, target_states[which]);

        return;

    };

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
        virus->get_prob_death(m) * (1.0 - p->get_death_reduction(virus, *m));

    // Recover
    m->array_double_tmp[1u] =
        1.0 - (1.0 - virus->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(virus, *m)); 
    

    // Running the roulette
    int which = roulette(2u, m);

    if (which < 0)
        return;

    // Which roulette happen?
    if (which == 0u) // If odd
    {

        int rm_state = -1;
        p->get_virus()->get_state(nullptr, nullptr, &rm_state);
        p->rm_virus(*m, rm_state);
        
    } else {

        p->rm_virus(*m);

    }

    return ;

};

#endif