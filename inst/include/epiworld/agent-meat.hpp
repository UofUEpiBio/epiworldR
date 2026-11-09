#ifndef EPIWORLD_PERSON_MEAT_HPP
#define EPIWORLD_PERSON_MEAT_HPP

#include <vector>
#include <string>
#include "config.hpp"
#include "epiworld-macros.hpp"
#include "agent-bones.hpp"
#include "agent-events-meat.hpp"

inline Agent::Agent(Agent && p) :
    state(p.state),
    state_prev(p.state_prev), 
    state_last_changed(p.state_last_changed),
    id(p.id),
    tools(std::move(p.tools)) /// Needs to be adjusted
{

    state = p.state;
    id     = p.id;
    
    // Dealing with the virus
    if (p.virus != nullptr)
    {
        virus = std::move(p.virus);
        virus->set_agent(this);
    }

    int loc = 0;
    for (auto & t : tools)
    {
        
        // Will create a copy of the virus, with the exeption of
        // the virus code
        t->agent     = this;
        t->pos_in_agent = loc++;

    }
    
}

// Copy constructor
inline Agent::Agent(const Agent & p)
{

    state              = p.state;
    state_prev         = p.state_prev;
    state_last_changed = p.state_last_changed;
    id                 = p.id;
    
    // Dealing with the virus
    if (p.virus != nullptr)
    {
        virus = std::shared_ptr<Virus>(p.virus->clone_ptr());
        virus->set_agent(this);
    }
    

    tools.clear();
    tools.reserve(p.get_n_tools());
    for (size_t i = 0u; i < p.get_n_tools(); ++i)
    {
        tools.emplace_back(std::shared_ptr<Tool>(p.tools[i]->clone_ptr()));
        tools.back()->set_agent(this, i);

    }
    
}

inline Agent & Agent::operator=(
    const Agent & other_agent
)
{

    state              = other_agent.state;
    state_prev         = other_agent.state_prev;
    state_last_changed = other_agent.state_last_changed;
    id                  = other_agent.id;
    
    if (other_agent.virus != nullptr)
    {
        virus = std::shared_ptr<Virus>(other_agent.virus->clone_ptr());
        virus->set_agent(this);
    } else
        virus = nullptr;
    
    
    tools.clear();
    tools.reserve(other_agent.get_n_tools());
    for (size_t i = 0u; i < other_agent.get_n_tools(); ++i)
    {
        tools.emplace_back(std::shared_ptr<Tool>(other_agent.tools[i]->clone_ptr()));
        tools.back()->set_agent(this, i);
    }
    
    return *this;
    
}

inline Agent::~Agent() = default;

inline void Agent::add_tool(Model & model, const Tool & tool)
{
    add_tool(model, tool, -99, -99);
}

inline void Agent::set_virus(Model & model, const Virus & virus)
{
    set_virus(model, virus, -99, -99);
}

inline void Agent::rm_tool(Model & model, epiworld_fast_uint tool_idx)
{
    rm_tool(model, tool_idx, -99, -99);
}

inline void Agent::rm_tool(Model & model, ToolPtr & tool)
{
    rm_tool(model, tool, -99, -99);
}

inline void Agent::add_tool(
    Model & model,
    const Tool & tool,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
) {

    // The model's own copy: only it has an id in the model (a tool given to
    // Model::add_tool() is copied, and the copy gets the id)
    if (
        (tool.get_id() < 0) ||
        (tool.get_id() >= static_cast<int>(model.get_db().get_n_tools()))
    )
        throw std::range_error(
            "The tool \"" + tool.get_name() + "\" is not in the model (id " +
            std::to_string(tool.get_id()) + "): give agents the model's copy, "
            "from Model::get_tool()."
        );

    ToolPtr tool_ptr = std::shared_ptr<Tool>(tool.clone_ptr());

    model._add_event(
        this, nullptr, tool_ptr, state_new, queue, EventAction::AddTool
        );

}

inline void Agent::set_virus(
    Model & model,
    const Virus & virus,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{
    
    // The model's own copy: only it has an id in the model (a virus given to
    // Model::add_virus() is copied, and the copy gets the id)
    if (
        (virus.get_id() < 0) ||
        (virus.get_id() >= static_cast<int>(model.get_db().get_n_viruses()))
    )
        throw std::range_error(
            "The virus \"" + virus.get_name() + "\" is not in the model (id " +
            std::to_string(virus.get_id()) + "): give agents the model's copy, "
            "from Model::get_virus()."
        );

    // With a natural history, a virus given from outside (seeding, an
    // import) follows the infection transition out of the agent's stage
    if (model.nh && (state_new == -99))
    {

        const int tr = model.nh->get_stage_info(model.nh_stage_of(state)).infection;
        if (tr < 0)
            throw std::logic_error(
                "Agent " + std::to_string(id) + " is in the stage \"" +
                model.nh->get_stages()[model.nh_stage_of(state)] + "\", which has no infection "
                "transition, so it cannot be given a virus."
            );

        state_new = model.nh->get_transition(tr).to;
        queue = Queue::NoOne;
        if (model.nh_pending == -1)
            model.nh_pending = Model::NH_SEED;

    }

    if (state_new == -99)
        throw std::logic_error(
            "A virus is given to agents of a model with a natural history."
        );

    if (queue == -99)
        queue = Queue::NoOne;

    VirusPtr virus_ptr = std::shared_ptr<Virus>(virus.clone_ptr());

    model._add_event(
        this, virus_ptr, nullptr, state_new, queue, EventAction::AddVirus,
        model.take_sampled_layer()
        );

}

inline void Agent::rm_tool(
    Model & model,
    epiworld_fast_uint tool_idx,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool_idx >= tools.size())
        throw std::range_error(
            "The Tool you want to remove is out of range. This Agent only has " +
            std::to_string(tools.size()) + " tools."
        );

    model._add_event(
        this, nullptr, tools[tool_idx], state_new, queue, EventAction::RemoveTool
        );

}

inline void Agent::rm_tool(
    Model & model,
    ToolPtr & tool,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (tool->agent != this)
        throw std::logic_error("Cannot remove a virus from another agent!");

    model._add_event(
        this, nullptr, tool, state_new, queue, EventAction::RemoveTool
        );

}

inline void Agent::rm_virus(
    Model & model,
    epiworld_fast_int state_new,
    epiworld_fast_int queue
)
{

    if (virus == nullptr)
        throw std::logic_error(
            "There is no virus to remove here!"
        );

    if (state_new == -99)
        throw std::logic_error(
            "Removing a virus needs the stage to move to: declare an external "
            "transition and fire() it."
        );

    if (queue == -99)
        queue = Queue::NoOne;

    model._add_event(
        this, virus, nullptr,
        state_new,
        queue,
        EventAction::RemoveVirus
        );

}

inline epiworld_double Agent::get_susceptibility_reduction(
    VirusPtr & v,
    Model & model
) {

    return model.susceptibility_reduction_mixer(this, v);
}

inline epiworld_double Agent::get_transmission_reduction(
    VirusPtr & v,
    Model & model
) {
    return model.transmission_reduction_mixer(this, v);
}

inline epiworld_double Agent::get_recovery_enhancer(
    VirusPtr & v,
    Model & model
) {
    return model.recovery_enhancer_mixer(this, v);
}

inline epiworld_double Agent::get_death_reduction(
    VirusPtr & v,
    Model & model
) {
    return model.death_reduction_mixer(this, v);
}

inline int Agent::get_id() const
{
    return id;
}

inline VirusPtr & Agent::get_virus() {
    return virus;
}

inline const VirusPtr & Agent::get_virus() const {
    return virus;
}


inline std::vector< ToolPtr > Agent::get_tools() {
    return tools;
}

inline const std::vector< ToolPtr > Agent::get_tools() const {
    return tools;
}

inline ToolPtr & Agent::get_tool(int i)
{
    return tools.at(i);
}

inline ToolPtr & Agent::get_tool(std::string name)
{
    for (auto & tool : tools)
        if (tool->get_name() == name)
            return tool;

    throw std::logic_error(
        "The agent does not have a tool with name: " + name
    );

}

inline size_t Agent::get_n_tools() const noexcept
{
    return tools.size();
}

inline void Agent::mutate_virus(Model & model)
{

    if (virus == nullptr)
        throw std::logic_error(
            "Agent " + std::to_string(id) + " has no virus to mutate."
        );

    virus->mutate(&model);

}

inline void Agent::change_state(
    Model & model,
    epiworld_fast_uint new_state,
    epiworld_fast_int queue
    )
{

    model._add_event(
        this, nullptr, nullptr, new_state, queue,
        EventAction::ChangeState
    );

    return;

}

inline unsigned int Agent::get_state() const {
    return state;
}

inline unsigned int Agent::get_state_prev() const {
    return state_prev;
}

inline int Agent::get_state_last_changed() const {
    return state_last_changed;
}

inline void Agent::reset()
{

    this->virus = nullptr;

    this->tools.clear();
    decltype(this->tools)().swap(this->tools);

    this->state = 0u;
    this->state_prev = 0u;

    this->state_last_changed = -1;
    
}

inline bool Agent::has_tool(epiworld_fast_uint t) const
{

    for (auto & tool : tools)
        if (tool->get_id() == static_cast<int>(t))
            return true;

    return false;

}

inline bool Agent::has_tool(std::string_view name) const
{

    for (auto & tool : tools)
        if (tool->get_name() == name)
            return true;

    return false;

}

inline bool Agent::has_tool(const Tool & tool) const
{

    return has_tool(tool.get_id());

}

inline bool Agent::has_virus(epiworld_fast_uint t) const
{
    if (virus->get_id() == static_cast<int>(t))
        return true;

    return false;
}

inline bool Agent::has_virus(std::string_view name) const
{
    
    if (virus->get_name() == name)
        return true;

    return false;

}

inline bool Agent::has_virus(const Virus & virus) const
{

    return has_virus(virus.get_id());

}

inline void Agent::print(
    Model & model,
    bool compressed
    ) const
{

    if (compressed)
    {
        printf_epiworld(
            "Agent: %i, state: %s (%i), Has virus: %s, NTools: %ii NNeigh: %i\n",
            static_cast<int>(id),
            model.states_labels[state].c_str(),
            static_cast<int>(state),
            virus == nullptr ? std::string("no").c_str() : std::string("yes").c_str(),
            static_cast<int>(tools.size()),
            static_cast<int>(model.network_degree(static_cast< size_t >(id)))
        );
    }
    else {

        printf_epiworld("Information about agent id %i\n",
            static_cast<int>(this->id));
        printf_epiworld("  State        : %s (%i)\n",
            model.states_labels[state].c_str(), static_cast<int>(state));
        printf_epiworld("  Has virus    : %s\n", virus == nullptr ?
            std::string("no").c_str() : std::string("yes").c_str());
        printf_epiworld("  Tool count   : %i\n", static_cast<int>(tools.size()));
        printf_epiworld("  Neigh. count : %i\n",
            static_cast<int>(model.network_degree(static_cast< size_t >(id))));

        size_t nfeats = model.get_agents_data_ncols();
        if (nfeats > 0)
        {

            printf_epiworld(
                "This model includes features (%i): [ ",
                static_cast<int>(nfeats)
                );

            int max_to_show = static_cast<int>((nfeats > 10)? 10 : nfeats);

            for (int k = 0; k < max_to_show; ++k)
            {
                printf_epiworld("%.2f", this->operator()(k, model));

                if (k != (max_to_show - 1))
                {
                    printf_epiworld(", ");
                } else {
                    printf_epiworld(" ]\n");
                }

            }
            
        }

    }

    return;

}

inline double & Agent::operator()(size_t j, Model & model)
{

    if (model.agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model.agents_data + j * model.size() + id);

}

inline double Agent::operator()(size_t j, const Model & model) const
{

    if (model.agents_data_ncols <= j)
        throw std::logic_error("The requested feature of the agent is out of range.");

    return *(model.agents_data + j * model.size() + id);

}

inline bool Agent::operator==(const Agent & other) const
{

    // Checking the address
    if (this == &other)
        return true;

    EPI_DEBUG_FAIL_AT_TRUE(
        state != other.state,
        "Agent:: state don't match"
        )
        

    EPI_DEBUG_FAIL_AT_TRUE(
        state_prev != other.state_prev,
        "Agent:: state_prev don't match"
        )
        

    // EPI_DEBUG_FAIL_AT_TRUE(
    //     state_last_changed != other.state_last_changed,
    //     "Agent:: state_last_changed don't match"
    //     ) ///< Last time the agent was updated.

    EPI_DEBUG_FAIL_AT_TRUE(
        ((virus == nullptr) && (other.virus != nullptr)) ||
            ((virus != nullptr) && (other.virus == nullptr)),
        "Agent:: virus don't match"
    )

    if ((virus != nullptr) && (other.virus != nullptr))
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *virus != *other.virus,
            "Agent:: virus doesn't match"
        )
    }
    
    EPI_DEBUG_FAIL_AT_TRUE(
        tools.size() != other.tools.size(),
        "Agent:: n_tools don't match"
    )

    for (size_t i = 0u; i < tools.size(); ++i)
    {
        
        EPI_DEBUG_FAIL_AT_TRUE(
            *tools[i] != *other.tools[i],
            "Agent:: tools[i] don't match"
        )
         
    }   
    
    return true;
    
}

#endif
