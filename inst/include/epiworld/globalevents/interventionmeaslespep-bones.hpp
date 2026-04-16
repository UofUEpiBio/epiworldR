#ifndef EPIWORLD_INTERVENTIONMEASLESPEP_BONES_HPP
#define EPIWORLD_INTERVENTIONMEASLESPEP_BONES_HPP

#include "../config.hpp"
#include "../model-bones.hpp"
#include "../agent-bones.hpp"

template<typename TSeq = EPI_DEFAULT_TSEQ>
class InterventionMeaslesPEP final : public GlobalEvent<TSeq> {

private:

    static inline const std::string _par_mmr_willingness{"PEP MMR willingness"};
    static inline const std::string _par_ig_willingness{"PEP IG willingness"};
    static inline const std::string _par_mmr_efficacy{"PEP MMR efficacy"};
    static inline const std::string _par_ig_efficacy{"PEP IG efficacy"};
    static inline const std::string _par_pep_mmr_window{"PEP MMR window"};
    static inline const std::string _par_pep_ig_window{"PEP IG window"};
    static inline const std::string _par_half_life_mean{"PEP IG half-life (mean)"};
    static inline const std::string _par_half_life_sd{"PEP IG half-life (sd)"};


    // Willigness and efficacy of the PEP
    epiworld_double _mmr_willingness;
    epiworld_double _ig_willingness;
    epiworld_double _mmr_efficacy;
    epiworld_double _ig_efficacy;
    epiworld_double _ig_half_life_mean;
    epiworld_double _ig_half_life_sd;
    epiworld_double _pep_mmr_window;
    epiworld_double _pep_ig_window;

    // Willingness of the agents to receive MMR and IG as PEP
    std::vector< bool > _willing_to_receive_mmr;
    std::vector< bool > _willing_to_receive_ig;

    // Target states to which the intervention applies
    std::vector< int > _target_states;
    std::vector< int > _states_if_pep_effective;
    std::vector< int > _states_if_pep_ineffective;

    // Id of the model
    int model_id = -1;

    /**
     * @brief Set up the intervention.
     * 
     * @details
     * This function is called at the beginning of the simulation. It sets up
     * the willingness of the agents to receive PEP and the efficacy of the
     * PEP.
     * @param model A pointer to the model.
     */
    void _setup(Model<TSeq> * model);

    // List of agents that will receive PEP
    // This is used to avoid iterating over the agents twice.
    std::vector< size_t > _to_receive_pep;
    std::vector< int > _next_if_effective;
    std::vector< int > _next_if_ineffective;

public:

    /**
     * @brief Construct a new Intervention PEP object.
     * @param par_mmr_efficacy The efficacy of the PEP. Must be between 0
     * and 1.
     * @param par_pep_willingness The willingness of the agents to receive
     * PEP. Must be between 0 and 1.
     * @param target_states The states to which the intervention applies. For
     * example, if the intervention applies to agents in quarantine, then this
     * should include the states that correspond to quarantine.
     * @param quarantine_states_for_pep The states to which the agents will be
     * moved if they receive PEP. Must be the same length as `quarantine_states`.
     * For example, if agents in state 2 are quarantined and will move to state
     * 5 if they receive PEP, then `quarantine_states` should include 2 and
     * `quarantine_states_for_pep` should include 5 at the corresponding
     * position.
     */
    InterventionMeaslesPEP(
        std::string name,
        epiworld_double mmr_efficacy,
        epiworld_double ig_efficacy,
        epiworld_double ig_half_life_mean,
        epiworld_double ig_half_life_sd,
        epiworld_double mmr_willingness,
        epiworld_double ig_willingness,
        epiworld_double mmr_window,
        epiworld_double ig_window,
        std::vector< int > target_states,
        std::vector< int > states_if_pep_effective,
        std::vector< int > states_if_pep_ineffective
    );

    /**
     * @brief Apply the intervention to the model.
     * @details
     * This function is called at the end of the day as a global event. It
     * iterates through the agents and gives PEP to those who are willing and
     * applicable.
     * 
     * Agents who receive PEP may then be moved to a different state if
     * the PEP is effective.
     */
    void operator()(Model<TSeq> * model, int day) override;

    std::unique_ptr< GlobalEvent<TSeq> > clone_ptr() const override;

    static bool agent_recovers(Agent<TSeq> & p, Model<TSeq> & m, int next_state);

};

#endif