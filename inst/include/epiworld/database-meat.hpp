#ifndef EPIWORLD_DATABASE_MEAT_HPP
#define EPIWORLD_DATABASE_MEAT_HPP

#include <vector>
#include <string>
#include "config.hpp"
#include "epiworld-macros.hpp"
#include "misc.hpp"
#include "database-bones.hpp"

/**
 * @file database-meat.hpp
 * @brief Recording: what the database does while a model runs. The
 * queries are in `database-meat-queries.hpp`, the output tables and files
 * in `database-meat-outputs.hpp`, and the comparison in
 * `database-meat-compare.hpp`.
 */

inline void DataBase::reset()
{

    // Initializing the counts
    rec.today_total.assign(model->nstates, 0);
    for (auto & p : model->get_agents())
        ++rec.today_total[p.get_state()];

    #ifdef EPI_DEBUG
    // Only the first should be different from zero
    {
        int n = static_cast<int>(model->size());
        if (rec.today_total[0] != n)
            throw std::runtime_error("The number of susceptible agents is not equal to the total number of agents.");

        if (std::accumulate(rec.today_total.begin(), rec.today_total.end(), 0) != n)
            throw std::runtime_error("The total number of agents is not equal to the sum of the number of agents in each state.");
            
    }
    #endif

    
    rec.transition_matrix.assign(model->nstates * model->nstates, 0);
    for (size_t s = 0u; s < model->nstates; ++s)
        rec.transition_matrix[s + s * model->nstates] = rec.today_total[s];

    // The last run's variants go; the founders keep their ids
    drop_variants();

    rec.hist_virus_date.clear();
    rec.hist_virus_id.clear();
    rec.hist_virus_state.clear();
    rec.hist_virus_counts.clear();

    rec.hist_tool_date.clear();
    rec.hist_tool_id.clear();
    rec.hist_tool_state.clear();
    rec.hist_tool_counts.clear();    

    rec.today_virus.resize(get_n_viruses());
    for (auto& virus_states : rec.today_virus)
        virus_states.assign(model->nstates, 0);

    rec.today_tool.resize(get_n_tools());
    for (auto& tool_states : rec.today_tool)
        tool_states.assign(model->nstates, 0);

    rec.hist_total_date.clear();
    rec.hist_total_state.clear();
    rec.hist_total_nviruses_active.clear();
    rec.hist_total_counts.clear();
    rec.hist_transition_matrix.clear();

    rec.transmission_date.clear();
    rec.transmission_virus.clear();
    rec.transmission_source.clear();
    rec.transmission_target.clear();
    rec.transmission_source_exposure_date.clear();
    rec.transmission_layer.clear();

    rec.transition_today.assign(
        !model->tracks.empty() ?
            model->tracks[0u].nh->get_n_transitions() +
                (model->st ? model->st->get_n_transitions() : 0u) :
            0u,
        0
    );
    rec.hist_transition_date.clear();
    rec.hist_transition_id.clear();
    rec.hist_transition_counts.clear();

    rec.m_hospitalizations.reset();

    return;

}

inline DataBase::DataBase(const DataBase & db) :
    reg(db.reg),
    rec(db.rec),
    sampling_freq(db.sampling_freq),
    user_data(nullptr)
{}

inline Model * DataBase::get_model() {
    return model;
}

inline void DataBase::record() 
{

    ////////////////////////////////////////////////////////////////////////////
    // DEBUGGING BLOCK
    ////////////////////////////////////////////////////////////////////////////
    EPI_DEBUG_SUM_INT(rec.today_total, model->size())
    EPI_DEBUG_ALL_NON_NEGATIVE(rec.today_total)

    #ifdef EPI_DEBUG
    // Checking whether the sums correspond
    std::vector< int > _today_total_cp(rec.today_total.size(), 0);
    for (auto & p : model->population)
        _today_total_cp[p.get_state()]++;
    
    EPI_DEBUG_VECTOR_MATCH_INT(
        _today_total_cp, rec.today_total,
        "Sums of __today_total_cp in database-meat.hpp"
        )

    // EPI_DEBUG only: the model's agents-by-state index (see
    // Model::get_agents_in_state()) is maintained incrementally, event by
    // event, so check it against the population once per step, next to the
    // same check for rec.today_total above. Each state's block must hold exactly
    // the agents in that state, each agent's recorded position must be right,
    // and the degree and carrier sums used by the transmission step must match
    // a recount. A drift here would silently bias which agents push, or which
    // mode "auto" picks.
    if (model->state_index_ready)
    {

        std::vector< size_t > _deg(model->nstates, 0u);
        std::vector< size_t > _carriers(model->nstates, 0u);
        std::vector< size_t > _carrier_deg(model->nstates, 0u);
        for (auto & p : model->population)
        {
            const size_t _d = model->network_degree(static_cast< size_t >(p.get_id()));
            _deg[p.get_state()] += _d;
            if (p.get_virus() != nullptr)
            {
                _carriers[p.get_state()]++;
                _carrier_deg[p.get_state()] += _d;
            }
        }

        for (size_t s = 0u; s < model->nstates; ++s)
        {

            const auto members = model->state_index_members(s);
            if (static_cast< int >(members.size()) != rec.today_total[s])
                throw std::logic_error("[epi-debug] DataBase::record state index size doesn't match rec.today_total.");

            for (size_t k = 0u; k < members.size(); ++k)
            {
                if (model->population[members[k]].get_state() != s)
                    throw std::logic_error("[epi-debug] DataBase::record state index lists an agent in another state.");
                if (model->agent_state[members[k]] != s)
                    throw std::logic_error("[epi-debug] DataBase::record agent_state is out of step.");
                if (model->agent_carrier[members[k]] != (model->population[members[k]].get_virus() != nullptr))
                    throw std::logic_error("[epi-debug] DataBase::record agent_carrier is out of step.");
                if (model->state_member_pos[members[k]] != model->state_start[s] + k)
                    throw std::logic_error("[epi-debug] DataBase::record state_member_pos is out of step.");
            }

            if (
                (_deg[s] != model->state_degree[s]) ||
                (_carriers[s] != model->state_carriers[s]) ||
                (_carrier_deg[s] != model->state_carrier_degree[s])
            )
                throw std::logic_error("[epi-debug] DataBase::record state degree/carrier sums are out of step.");

        }

    }

    if (model->today() == 0)
    {
        if (rec.hist_total_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_total_date should be of length 0.")
        if (rec.hist_total_nviruses_active.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_total_nviruses_active should be of length 0.")
        if (rec.hist_total_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_total_state should be of length 0.")
        if (rec.hist_total_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_total_counts should be of length 0.")
        if (rec.hist_virus_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_virus_date should be of length 0.")
        if (rec.hist_virus_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_virus_id should be of length 0.")
        if (rec.hist_virus_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_virus_state should be of length 0.")
        if (rec.hist_virus_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_virus_counts should be of length 0.")
        if (rec.hist_tool_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_tool_date should be of length 0.")
        if (rec.hist_tool_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_tool_id should be of length 0.")
        if (rec.hist_tool_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_tool_state should be of length 0.")
        if (rec.hist_tool_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record rec.hist_tool_counts should be of length 0.")
    }
    #endif
    ////////////////////////////////////////////////////////////////////////////

    // Only store every now and then
    if ((model->today() % sampling_freq) == 0)
    {

        // Recording virus's history
        for (size_t v = 0u; v < rec.today_virus.size(); ++v)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                rec.hist_virus_date.push_back(model->today());
                rec.hist_virus_id.push_back(static_cast< int >(v));
                rec.hist_virus_state.push_back(s);

                #ifdef EPI_DEBUG
                if (rec.today_virus[v][s] < 0)
                    throw std::logic_error(
                        "The count of viruses in DataBase::record cannot be negative."
                    );
                #endif

                rec.hist_virus_counts.push_back(rec.today_virus[v][s]);

            }

        }

        // Recording tool's history
        for (size_t t = 0u; t < rec.today_tool.size(); ++t)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                rec.hist_tool_date.push_back(model->today());
                rec.hist_tool_id.push_back(static_cast< int >(t));
                rec.hist_tool_state.push_back(s);
                rec.hist_tool_counts.push_back(rec.today_tool[t][s]);

            }

        }

        // Recording the overall history
        for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
        {
            rec.hist_total_date.push_back(model->today());
            rec.hist_total_nviruses_active.push_back(rec.today_total_nviruses_active);
            rec.hist_total_state.push_back(s);
            rec.hist_total_counts.push_back(rec.today_total[s]);
        }

        for (const auto& cell : rec.transition_matrix)
            rec.hist_transition_matrix.push_back(cell);

        for (size_t k = 0u; k < rec.transition_today.size(); ++k)
        {
            rec.hist_transition_date.push_back(model->today());
            rec.hist_transition_id.push_back(static_cast< int >(k));
            rec.hist_transition_counts.push_back(rec.transition_today[k]);
            rec.transition_today[k] = 0;
        }

        // Now the diagonal must reflect the state
        for (size_t s_i = 0u; s_i < model->nstates; ++s_i)
        {

            for (size_t s_j = 0u; s_j < model->nstates; ++s_j)
            {
                
                if ((s_i != s_j) && (rec.transition_matrix[s_i + s_j * model->nstates] > 0))
                {
                    rec.transition_matrix[s_j + s_j * model->nstates] +=
                        rec.transition_matrix[s_i + s_j * model->nstates];

                    rec.transition_matrix[s_i + s_j * model->nstates] = 0;
                }
         
            }

        }

        #ifdef EPI_DEBUG
        for (size_t s_i = 0u; s_i < model->nstates; ++s_i)
        {
            if (rec.transition_matrix[s_i + s_i * model->nstates] != 
                rec.today_total[s_i])
                throw std::logic_error(
                    "The diagonal of the updated transition Matrix should match the daily totals"
                    );
        }
        #endif

    }

}

inline void DataBase::drop_variants()
{

    const size_t n = reg.n_founders;
    if (reg.virus_name.size() <= n)
        return;

    reg.virus_name.resize(n);
    reg.virus_origin_date.resize(n);
    reg.virus_parent_id.resize(n);
    rec.today_virus.resize(n);
    rec.today_total_nviruses_active = static_cast< int >(n);

}

inline void DataBase::record_virus(Virus & v, bool founder)
{

    if (founder)
    {
        drop_variants();
        ++reg.n_founders;
    }
    else if (v.get_id() < 0)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" mutated but is not on record."
        );

    const int parent = founder ? -99 : v.get_id();
    const int id = static_cast< int >(reg.virus_name.size());

    reg.virus_name.push_back(v.get_name());
    reg.virus_origin_date.push_back(model->today());
    reg.virus_parent_id.push_back(parent);
    rec.today_virus.push_back(std::vector< int >(model->nstates, 0));
    rec.today_total_nviruses_active++;

    v.set_id(id);
    v.set_date(model->today());

    // A founder starts a lineage, which its mutations keep. Lineages 63 and
    // above share the overflow bit, which tools cannot target.
    if (founder)
    {
        v.lineage_id  = id;
        v.lineage_bit = uint64_t(1) << std::min(id, 63);
    }
    else if (v.get_agent() != nullptr)
    {
        // The agent's counts move from the parent to the variant
        const auto s = v.get_agent()->get_state();
        rec.today_virus[parent][s]--;
        rec.today_virus[id][s]++;
    }

}

inline void DataBase::record_tool(Tool & t)
{

    t.set_id(static_cast< int >(reg.tool_name.size()));
    t.set_date(model->today());
    reg.tool_name.push_back(t.get_name());
    reg.tool_origin_date.push_back(model->today());
    rec.today_tool.push_back(std::vector< int >(model->nstates, 0));

}

inline size_t DataBase::size() const
{
    return reg.virus_name.size();
}

inline void DataBase::update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    if (prev_state == new_state)
        return; // No need to update if the state is the same

    #ifdef EPI_DEBUG
    // Checking ranges (should be within expected)
    if ((prev_state >= model->nstates))
        throw std::out_of_range(
            "prev_state is out of range in DataBase::update_state"
        );
    if ((new_state >= model->nstates))
        throw std::out_of_range(
            "new_state is out of range in DataBase::update_state"
        );
    if (prev_state == new_state)
        throw std::logic_error(
            "prev_state and new_state are the same in DataBase::update_state"
        );
    #endif

    if (undo)
    {

        rec.today_total[prev_state]++;
        rec.today_total[new_state]--;
        
    } else {

        rec.today_total[prev_state]--;
        rec.today_total[new_state]++;

    }

    record_transition(prev_state, new_state, undo);
    
    return;
}

inline void DataBase::update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    if (undo)
    {
        rec.today_virus[virus_id][prev_state]++;
        rec.today_virus[virus_id][new_state]--;
    }
    else
    {
        rec.today_virus[virus_id][prev_state]--;
        rec.today_virus[virus_id][new_state]++;
    }

    #ifdef EPI_DEBUG
    if (rec.today_virus[virus_id][prev_state] < 0)
        throw std::logic_error("An entry in rec.today_virus is negative.");
    if (rec.today_virus[virus_id][new_state] < 0)
        throw std::logic_error("An entry in rec.today_virus is negative (new_state).");
    #endif

    return;
    
}

inline void DataBase::update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    if (undo)
    {
        rec.today_tool[tool_id][prev_state]++;
        rec.today_tool[tool_id][new_state]--;
    }
    else
    {
        rec.today_tool[tool_id][prev_state]--;
        rec.today_tool[tool_id][new_state]++;
    }

    return;

}

inline void DataBase::record_transition(
    epiworld_fast_uint from,
    epiworld_fast_uint to,
    bool undo
) {

    if (undo)
    {   

        rec.transition_matrix[to * model->nstates + from]--;
        rec.transition_matrix[from * model->nstates + from]++;

    } else {

        rec.transition_matrix[to * model->nstates + from]++;
        rec.transition_matrix[from * model->nstates + from]--;

    }

    #ifdef EPI_DEBUG
    if (rec.transition_matrix[from * model->nstates + from] < 0)
        throw std::logic_error("An entry in transition matrix is negative.");
    #endif

}

inline void DataBase::record_transmission(
    int i,
    int j,
    int virus,
    int i_expo_date,
    int layer
) {

    rec.transmission_date.push_back(model->today());
    rec.transmission_source.push_back(i);
    rec.transmission_target.push_back(j);
    rec.transmission_virus.push_back(virus);
    rec.transmission_source_exposure_date.push_back(i_expo_date);
    rec.transmission_layer.push_back(layer);

}

inline size_t DataBase::get_n_viruses() const
{
    return reg.virus_name.size();
}

inline size_t DataBase::get_n_tools() const
{
    return reg.tool_name.size();
}

inline void DataBase::set_user_data(
    std::vector< std::string > names
)
{
    user_data = UserData(names);
    user_data.model = model;
}

inline void DataBase::add_user_data(
    std::vector< epiworld_double > x
)
{

    user_data.add(x);

}

inline void DataBase::add_user_data(
    epiworld_fast_uint k,
    epiworld_double x
)
{

    user_data.add(k, x);

}

inline UserData & DataBase::get_user_data()
{
    return user_data;
}

inline void DataBase::record_hospitalization(Agent & agent)
{
    rec.m_hospitalizations.record(agent, *model);
}

#endif
