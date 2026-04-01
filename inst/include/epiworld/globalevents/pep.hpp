#ifndef EPIWORLD_MODELS_INTERVENTIONS_HPP
#define EPIWORLD_MODELS_INTERVENTIONS_HPP

#include "../config.hpp"
#include "../model-bones.hpp"
#include "../agent-bones.hpp"
#include "../tools/vaccine.hpp"

template<typename TSeq = EPI_DEFAULT_TSEQ>
class InterventionPEP : public GlobalEvent<TSeq> {

private: 

    // Willigness and efficacy of the PEP
    std::string _parname_willingness;
    std::string _parname_efficacy;
    
    // Willingness of the agents to receive PEP
    std::vector< bool > _willing_to_receive_pep;

    // Quarantine states to which the intervention applies
    std::vector< int > _quarantine_states;
    std::vector< int > _quarantine_states_for_pep;

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

public:

    /**
     * @brief Construct a new Intervention PEP object.
     * @param parameter_efficacy The efficacy of the PEP. Must be between 0
     * and 1.
     * @param parameter_willingness The willingness of the agents to receive
     * PEP. Must be between 0 and 1.
     * @param quarantine_states The states to which the intervention applies. For
     * example, if the intervention applies to agents in quarantine, then this
     * should include the states that correspond to quarantine.
     */
    InterventionPEP(
        std::string parameter_efficacy,
        std::string parameter_willingness,
        std::vector< int > quarantine_states,
        std::vector< int > quarantine_states_for_pep
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

template<typename TSeq>
inline InterventionPEP<TSeq>::InterventionPEP(
    std::string parameter_efficacy,
    std::string parameter_willingness,
    std::vector< int > quarantine_states,
    std::vector< int > quarantine_states_for_pep
) {

    // Must match the length
    if (quarantine_states.size() != quarantine_states_for_pep.size())
        throw std::logic_error(
            "The length of the quarantine states and the quarantine states for "
            "PEP must be the same. These are currently: " +
            std::to_string(quarantine_states.size()) + " and " +
            std::to_string(quarantine_states_for_pep.size()) +
            ", respectively."
        );

    this->_quarantine_states = quarantine_states;
    this->_quarantine_states_for_pep = quarantine_states_for_pep;
    this->_parname_efficacy = parameter_efficacy;
    this->_parname_willingness = parameter_willingness;
}

template<typename TSeq>
inline void InterventionPEP<TSeq>::_setup(
    Model<TSeq> * model
) {

    // Randomizing willingness
    this->_willing_to_receive_pep.assign(model->size(), false);

    auto willingness = model->par(this->_parname_willingness);
    for (size_t i = 0u; i < model->size(); ++i)
    {
        if (model->runif() < willingness)
        {
            this->_willing_to_receive_pep[i] = true;
        }
    }

    // Adding the PEP vaccine as a tool to the model
    // (if not already added by the user)
    if (model->has_tool("PEP MMR"))
        return;

    // Creating the PEP vaccine tool
    ToolVaccine<TSeq> pep("PEP MMR");
    pep.set_susceptibility_reduction(model->par(this->_parname_efficacy));
    model->add_tool(pep);

};

template<typename TSeq>
inline void InterventionPEP<TSeq>::operator()(Model<TSeq> * model, int) {

    // Verifying if this needs to be setup
    if (static_cast<int>(model->get_sim_id()) != this->model_id)
    {
        this->model_id = static_cast<int>(model->get_sim_id());
        this->_setup(model);
    }

    // Iterate through the agents and check 
    for (auto & agent: model->get_agents())
    {

        // Checking if the agent is in a quarantine
        // state
        int agent_state = static_cast<int>(agent.get_state());
        if (!IN(agent_state, this->_quarantine_states))
            continue;

        // Checking willigness
        if (!this->_willing_to_receive_pep[agent.get_id()])
            continue;

        // Finding the corresponding state for PEP
        auto it = std::find(
            this->_quarantine_states.begin(),
            this->_quarantine_states.end(),
            agent.get_state()
        );

        // No need to check it, we know it is there
        auto pos = std::distance(this->_quarantine_states.begin(), it);

        // We will administer PEP to the agent
        agent.add_tool(
            *model,
            model->get_tool("PEP MMR"),
            this->_quarantine_states_for_pep[pos]
        );

    }

};

template<typename TSeq>
inline std::unique_ptr< GlobalEvent<TSeq> > InterventionPEP<TSeq>::clone_ptr() const
{
    return std::make_unique< InterventionPEP<TSeq>>(*this);

}

template<typename TSeq>
inline bool InterventionPEP<TSeq>::agent_recovers(
    Agent<TSeq> & p,
    Model<TSeq> & m,
    int next_state
) {

    // If the agent has PEP, then we have to figure out if it works or not
    if (p.has_tool("PEP MMR"))
    {

        // Adding the tool
        auto & tool = p.get_tool("PEP MMR");
        if (tool.get_susceptibility_reduction(p.get_virus(), m) > 0.0)
        {
            p.rm_virus(m, next_state);
            return true;
        }
    }

    return false;

}

#endif