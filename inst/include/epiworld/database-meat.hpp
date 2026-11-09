#ifndef EPIWORLD_DATABASE_MEAT_HPP
#define EPIWORLD_DATABASE_MEAT_HPP

#include <vector>
#include <string>
#include "config.hpp"
#include "epiworld-macros.hpp"
#include "misc.hpp"
#include "database-bones.hpp"

inline void DataBase::reset()
{

    // Initializing the counts
    today_total.assign(model->nstates, 0);
    for (auto & p : model->get_agents())
        ++today_total[p.get_state()];

    #ifdef EPI_DEBUG
    // Only the first should be different from zero
    {
        int n = static_cast<int>(model->size());
        if (today_total[0] != n)
            throw std::runtime_error("The number of susceptible agents is not equal to the total number of agents.");

        if (std::accumulate(today_total.begin(), today_total.end(), 0) != n)
            throw std::runtime_error("The total number of agents is not equal to the sum of the number of agents in each state.");
            
    }
    #endif

    
    transition_matrix.assign(model->nstates * model->nstates, 0);
    for (size_t s = 0u; s < model->nstates; ++s)
        transition_matrix[s + s * model->nstates] = today_total[s];

    // The last run's variants go; the founders keep their ids
    drop_variants();

    hist_virus_date.clear();
    hist_virus_id.clear();
    hist_virus_state.clear();
    hist_virus_counts.clear();

    hist_tool_date.clear();
    hist_tool_id.clear();
    hist_tool_state.clear();
    hist_tool_counts.clear();    

    today_virus.resize(get_n_viruses());
    for (auto& virus_states : today_virus)
        virus_states.assign(model->nstates, 0);

    today_tool.resize(get_n_tools());
    for (auto& tool_states : today_tool)
        tool_states.assign(model->nstates, 0);

    hist_total_date.clear();
    hist_total_state.clear();
    hist_total_nviruses_active.clear();
    hist_total_counts.clear();
    hist_transition_matrix.clear();

    transmission_date.clear();
    transmission_virus.clear();
    transmission_source.clear();
    transmission_target.clear();
    transmission_source_exposure_date.clear();
    transmission_layer.clear();

    transition_today.assign(
        model->nh ?
            model->nh->get_n_transitions() +
                (model->st ? model->st->get_n_transitions() : 0u) :
            0u,
        0
    );
    hist_transition_date.clear();
    hist_transition_id.clear();
    hist_transition_counts.clear();

    m_hospitalizations.reset();

    return;

}

inline DataBase::DataBase(const DataBase & db) :
    m_hospitalizations(db.m_hospitalizations),
    n_founders(db.n_founders),
    virus_name(db.virus_name),
    virus_origin_date(db.virus_origin_date),
    virus_parent_id(db.virus_parent_id),
    tool_name(db.tool_name),
    tool_origin_date(db.tool_origin_date),
    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    today_virus(db.today_virus),
    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    today_tool(db.today_tool),
    // {Susceptible, Infected, etc.}
    today_total(db.today_total),
    // Totals
    today_total_nviruses_active(db.today_total_nviruses_active),
    sampling_freq(db.sampling_freq),
    // Variants history
    hist_virus_date(db.hist_virus_date),
    hist_virus_id(db.hist_virus_id),
    hist_virus_state(db.hist_virus_state),
    hist_virus_counts(db.hist_virus_counts),
    // Tools history
    hist_tool_date(db.hist_tool_date),
    hist_tool_id(db.hist_tool_id),
    hist_tool_state(db.hist_tool_state),
    hist_tool_counts(db.hist_tool_counts),
    // Overall hist
    hist_total_date(db.hist_total_date),
    hist_total_nviruses_active(db.hist_total_nviruses_active),
    hist_total_state(db.hist_total_state),
    hist_total_counts(db.hist_total_counts),
    hist_transition_matrix(db.hist_transition_matrix),
    // Transmission network
    transmission_date(db.transmission_date),
    transmission_source(db.transmission_source),
    transmission_target(db.transmission_target),
    transmission_virus(db.transmission_virus),
    transmission_source_exposure_date(db.transmission_source_exposure_date),
    transmission_layer(db.transmission_layer),
    transition_matrix(db.transition_matrix),
    transition_today(db.transition_today),
    hist_transition_date(db.hist_transition_date),
    hist_transition_id(db.hist_transition_id),
    hist_transition_counts(db.hist_transition_counts),
    user_data(nullptr)
{}

// DataBase & DataBase::operator=(const DataBase & m)
// {

// }

inline Model * DataBase::get_model() {
    return model;
}

inline void DataBase::record() 
{

    ////////////////////////////////////////////////////////////////////////////
    // DEBUGGING BLOCK
    ////////////////////////////////////////////////////////////////////////////
    EPI_DEBUG_SUM_INT(today_total, model->size())
    EPI_DEBUG_ALL_NON_NEGATIVE(today_total)

    #ifdef EPI_DEBUG
    // Checking whether the sums correspond
    std::vector< int > _today_total_cp(today_total.size(), 0);
    for (auto & p : model->population)
        _today_total_cp[p.get_state()]++;
    
    EPI_DEBUG_VECTOR_MATCH_INT(
        _today_total_cp, today_total,
        "Sums of __today_total_cp in database-meat.hpp"
        )

    // EPI_DEBUG only: the model's agents-by-state index (see
    // Model::get_agents_in_state()) is maintained incrementally, event by
    // event, so check it against the population once per step, next to the
    // same check for today_total above. Each state's block must hold exactly
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
            if (static_cast< int >(members.size()) != today_total[s])
                throw std::logic_error("[epi-debug] DataBase::record state index size doesn't match today_total.");

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
        if (hist_total_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_date should be of length 0.")
        if (hist_total_nviruses_active.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_nviruses_active should be of length 0.")
        if (hist_total_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_state should be of length 0.")
        if (hist_total_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_total_counts should be of length 0.")
        if (hist_virus_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_virus_date should be of length 0.")
        if (hist_virus_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_virus_id should be of length 0.")
        if (hist_virus_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_virus_state should be of length 0.")
        if (hist_virus_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_virus_counts should be of length 0.")
        if (hist_tool_date.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_date should be of length 0.")
        if (hist_tool_id.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_id should be of length 0.")
        if (hist_tool_state.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_state should be of length 0.")
        if (hist_tool_counts.size() != 0)
            EPI_DEBUG_ERROR(std::logic_error, "DataBase::record hist_tool_counts should be of length 0.")
    }
    #endif
    ////////////////////////////////////////////////////////////////////////////

    // Only store every now and then
    if ((model->today() % sampling_freq) == 0)
    {

        // Recording virus's history
        for (size_t v = 0u; v < today_virus.size(); ++v)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                hist_virus_date.push_back(model->today());
                hist_virus_id.push_back(static_cast< int >(v));
                hist_virus_state.push_back(s);

                #ifdef EPI_DEBUG
                if (today_virus[v][s] < 0)
                    throw std::logic_error(
                        "The count of viruses in DataBase::record cannot be negative."
                    );
                #endif

                hist_virus_counts.push_back(today_virus[v][s]);

            }

        }

        // Recording tool's history
        for (size_t t = 0u; t < today_tool.size(); ++t)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                hist_tool_date.push_back(model->today());
                hist_tool_id.push_back(static_cast< int >(t));
                hist_tool_state.push_back(s);
                hist_tool_counts.push_back(today_tool[t][s]);

            }

        }

        // Recording the overall history
        for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
        {
            hist_total_date.push_back(model->today());
            hist_total_nviruses_active.push_back(today_total_nviruses_active);
            hist_total_state.push_back(s);
            hist_total_counts.push_back(today_total[s]);
        }

        for (const auto& cell : transition_matrix)
            hist_transition_matrix.push_back(cell);

        for (size_t k = 0u; k < transition_today.size(); ++k)
        {
            hist_transition_date.push_back(model->today());
            hist_transition_id.push_back(static_cast< int >(k));
            hist_transition_counts.push_back(transition_today[k]);
            transition_today[k] = 0;
        }

        // Now the diagonal must reflect the state
        for (size_t s_i = 0u; s_i < model->nstates; ++s_i)
        {

            for (size_t s_j = 0u; s_j < model->nstates; ++s_j)
            {
                
                if ((s_i != s_j) && (transition_matrix[s_i + s_j * model->nstates] > 0))
                {
                    transition_matrix[s_j + s_j * model->nstates] +=
                        transition_matrix[s_i + s_j * model->nstates];

                    transition_matrix[s_i + s_j * model->nstates] = 0;
                }
         
            }

        }

        #ifdef EPI_DEBUG
        for (size_t s_i = 0u; s_i < model->nstates; ++s_i)
        {
            if (transition_matrix[s_i + s_i * model->nstates] != 
                today_total[s_i])
                throw std::logic_error(
                    "The diagonal of the updated transition Matrix should match the daily totals"
                    );
        }
        #endif

    }

}

inline void DataBase::drop_variants()
{

    const size_t n = n_founders;
    if (virus_name.size() <= n)
        return;

    virus_name.resize(n);
    virus_origin_date.resize(n);
    virus_parent_id.resize(n);
    today_virus.resize(n);
    today_total_nviruses_active = static_cast< int >(n);

}

inline void DataBase::record_virus(Virus & v, bool founder)
{

    if (founder)
    {
        drop_variants();
        ++n_founders;
    }
    else if (v.get_id() < 0)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" mutated but is not on record."
        );

    const int parent = founder ? -99 : v.get_id();
    const int id = static_cast< int >(virus_name.size());

    virus_name.push_back(v.get_name());
    virus_origin_date.push_back(model->today());
    virus_parent_id.push_back(parent);
    today_virus.push_back(std::vector< int >(model->nstates, 0));
    today_total_nviruses_active++;

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
        today_virus[parent][s]--;
        today_virus[id][s]++;
    }

}

inline void DataBase::record_tool(Tool & t)
{

    t.set_id(static_cast< int >(tool_name.size()));
    t.set_date(model->today());
    tool_name.push_back(t.get_name());
    tool_origin_date.push_back(model->today());
    today_tool.push_back(std::vector< int >(model->nstates, 0));

}

inline size_t DataBase::size() const
{
    return virus_name.size();
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

        today_total[prev_state]++;
        today_total[new_state]--;
        
    } else {

        today_total[prev_state]--;
        today_total[new_state]++;

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
        today_virus[virus_id][prev_state]++;
        today_virus[virus_id][new_state]--;
    }
    else
    {
        today_virus[virus_id][prev_state]--;
        today_virus[virus_id][new_state]++;
    }

    #ifdef EPI_DEBUG
    if (today_virus[virus_id][prev_state] < 0)
        throw std::logic_error("An entry in today_virus is negative.");
    if (today_virus[virus_id][new_state] < 0)
        throw std::logic_error("An entry in today_virus is negative (new_state).");
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
        today_tool[tool_id][prev_state]++;
        today_tool[tool_id][new_state]--;
    }
    else
    {
        today_tool[tool_id][prev_state]--;
        today_tool[tool_id][new_state]++;
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

        transition_matrix[to * model->nstates + from]--;
        transition_matrix[from * model->nstates + from]++;

    } else {

        transition_matrix[to * model->nstates + from]++;
        transition_matrix[from * model->nstates + from]--;

    }

    #ifdef EPI_DEBUG
    if (transition_matrix[from * model->nstates + from] < 0)
        throw std::logic_error("An entry in transition matrix is negative.");
    #endif

}

inline void DataBase::split_state(
    epiworld_fast_uint state,
    std::string & status,
    std::string & stage
) const
{

    if (!model->nh)
    {
        status = "Free";
        stage = model->states_labels[state];
        return;
    }

    const unsigned int c = static_cast< unsigned int >(state);
    stage = model->nh->get_stages()[static_cast< size_t >(model->nh_stage_of(c))];
    status = model->st ?
        model->st->get_stages()[static_cast< size_t >(model->nh_status_of(c))] :
        std::string("Free");

}

inline void DataBase::split_states(
    const std::vector< epiworld_fast_uint > & states,
    std::vector< std::string > * status,
    std::vector< std::string > * stage
) const
{

    std::string k, s;
    if (status != nullptr)
        status->resize(states.size());
    if (stage != nullptr)
        stage->resize(states.size());

    for (size_t i = 0u; i < states.size(); ++i)
    {
        split_state(states[i], k, s);
        if (status != nullptr)
            (*status)[i] = k;
        if (stage != nullptr)
            (*stage)[i] = s;
    }

}

inline int DataBase::today_total_of(
    const std::string & name,
    bool by_stage
) const
{

    bool found = false;
    int total = 0;
    std::string k, s;
    for (epiworld_fast_uint i = 0u; i < today_total.size(); ++i)
    {
        split_state(i, k, s);
        if ((by_stage ? s : k) == name)
        {
            found = true;
            total += today_total[i];
        }
    }

    if (!found)
        throw std::range_error(
            "'" + name + "' is not a " + (by_stage ? "stage" : "status") +
            " of the model."
        );

    return total;

}

inline int DataBase::get_today_total_stage(
    const std::string & stage
) const
{
    return today_total_of(stage, true);
}

inline int DataBase::get_today_total_status(
    const std::string & status
) const
{
    return today_total_of(status, false);
}

inline void DataBase::get_today_total(
    std::vector< std::string > * status,
    std::vector< std::string > * stage,
    std::vector< int > * counts
) const
{

    std::vector< epiworld_fast_uint > states(today_total.size());
    for (size_t i = 0u; i < states.size(); ++i)
        states[i] = static_cast< epiworld_fast_uint >(i);

    split_states(states, status, stage);

    if (counts != nullptr)
        *counts = today_total;

}

inline void DataBase::get_today_virus(
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & id,
    std::vector< int > & counts
    ) const
{

    // One row per virus and state
    const size_t n_viruses = today_virus.size();
    const size_t n_states = model->states_labels.size();

    std::vector< epiworld_fast_uint > states;
    states.reserve(n_viruses * n_states);
    id.clear();
    counts.clear();

    for (epiworld_fast_uint v = 0u; v < n_viruses; ++v)
        for (epiworld_fast_uint s = 0u; s < n_states; ++s)
        {
            states.push_back(s);
            id.push_back(static_cast<int>(v));
            counts.push_back(today_virus[v][s]);
        }

    split_states(states, &status, &stage);

}

inline void DataBase::get_hist_total(
    std::vector< int > * date,
    std::vector< std::string > * status,
    std::vector< std::string > * stage,
    std::vector< int > * counts
) const
{

    if (date != nullptr)
        *date = hist_total_date;

    split_states(hist_total_state, status, stage);

    if (counts != nullptr)
        *counts = hist_total_counts;

}

inline void DataBase::get_hist_virus(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & counts
) const {

    date = hist_virus_date;
    id = hist_virus_id;
    split_states(hist_virus_state, &status, &stage);
    counts = hist_virus_counts;

}

inline void DataBase::get_hist_tool(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & counts
) const {

    date = hist_tool_date;
    id = hist_tool_id;
    split_states(hist_tool_state, &status, &stage);
    counts = hist_tool_counts;

}

inline void DataBase::get_today_transition_matrix(
    std::vector< int > & counts
) const
{

    counts = transition_matrix;

    return;

}

inline void DataBase::get_hist_transition_matrix(
    std::vector< std::string > & status_from,
    std::vector< std::string > & stage_from,
    std::vector< std::string > & status_to,
    std::vector< std::string > & stage_to,
    std::vector< int > & date,
    std::vector< int > & counts,
    bool skip_zeros
) const
{

    size_t n = this->hist_transition_matrix.size();

    std::vector< epiworld_fast_uint > from, to;
    from.reserve(n);
    to.reserve(n);
    date.clear();
    counts.clear();
    date.reserve(n);
    counts.reserve(n);

    size_t n_states = model->nstates;
    size_t n_steps  = model->get_ndays();

    // If n is zero, then we are done
    if (n == 0u)
    {
        split_states(from, &status_from, &stage_from);
        split_states(to, &status_to, &stage_to);
        return;
    }

    for (size_t step = 0u; step <= n_steps; ++step) // The final step counts
    {
        for (size_t j = 0u; j < n_states; ++j) // Column major storage
        {
            for (size_t i = 0u; i < n_states; ++i)
            {
                // Retrieving the value of the day
                int v = hist_transition_matrix[
                    step * n_states * n_states + // Day of the data
                    j * n_states +               // Column (to)
                    i                            // Row (from)
                    ];

                // If we are skipping the zeros and it is zero, then don't save
                if (skip_zeros && v == 0)
                    continue;
                                
                from.push_back(static_cast< epiworld_fast_uint >(i));
                to.push_back(static_cast< epiworld_fast_uint >(j));
                date.push_back(hist_total_date[step * n_states]);
                counts.push_back(v);

            }

        }
    }

    split_states(from, &status_from, &stage_from);
    split_states(to, &status_to, &stage_to);


}


inline void DataBase::get_active_cases(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & count
) const 
{

    // Extracting useful sizes
    size_t n_days    = model->get_ndays() + 1u; // Including day 0
    size_t n_viruses = model->get_n_viruses();
    
    // Making room
    date.assign(n_days * n_viruses, 0);
    virus_id.assign(n_days * n_viruses, 0);
    count.assign(n_days * n_viruses, 0);

    for (size_t i = 0u; i < hist_virus_id.size(); ++i)
    {

        // With more viruses, there are more days
        auto location = hist_virus_date[i] + hist_virus_id[i] * n_days;

        date[location] = hist_virus_date[i];
        virus_id[location] = hist_virus_id[i];
        count[location] += hist_virus_counts[i];

    }

    return;
    
}

inline void DataBase::get_outbreak_size(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & outbreak_size
) const 
{

    // Extracting useful sizes
    size_t n_days    = model->get_ndays() + 1u; // Including day 0
    size_t n_viruses = model->get_n_viruses();
    
    // Making room
    date.assign(n_days * n_viruses, 0);
    std::iota(date.begin(), date.end(), 0);

    virus_id.assign(n_days * n_viruses, 0);
    for (size_t v = 0u; v < n_viruses; ++v)
        for (size_t d = 0u; d < n_days; ++d)
            virus_id[d + v * n_days] = static_cast<int>(v);

    outbreak_size.assign(n_days * n_viruses, 0);

    // With more viruses, there are more days;
    for (size_t i = 0u; i < transmission_date.size(); ++i)
    {
        outbreak_size[
            transmission_date[i] + transmission_virus[i] * n_days
        ] += 1;
    }

    // Now, we generate the cumulative sum
    for (size_t v = 0u; v < n_viruses; ++v)
    {
        for (size_t d = 1u; d < n_days; ++d)
        {
            auto location = d + v * n_days;
            outbreak_size[location] += outbreak_size[location - 1u];
        }
    }

    return;
    
}

inline void DataBase::get_transmissions(
    std::vector<int> & date,
    std::vector<int> & source,
    std::vector<int> & target,
    std::vector<int> & virus,
    std::vector<int> & source_exposure_date,
    std::vector<int> & layer
) const 
{

    size_t nevents = transmission_date.size();

    date.resize(nevents);
    source.resize(nevents);
    target.resize(nevents);
    virus.resize(nevents);
    source_exposure_date.resize(nevents);
    layer.resize(nevents);

    get_transmissions(
        date.data(),
        source.data(),
        target.data(),
        virus.data(),
        source_exposure_date.data(),
        layer.data()
    );

}

inline void DataBase::get_transmissions(
    int * date,
    int * source,
    int * target,
    int * virus,
    int * source_exposure_date,
    int * layer
) const 
{

    size_t nevents = transmission_date.size();

    for (size_t i = 0u; i < nevents; ++i)
    {

        *(date + i) = transmission_date.at(i);
        *(source + i) = transmission_source.at(i);
        *(target + i) = transmission_target.at(i);
        *(virus + i) = transmission_virus.at(i);
        *(source_exposure_date + i) = transmission_source_exposure_date.at(i);
        *(layer + i) = transmission_layer.at(i);

    }

}

inline RunOutputs DataBase::get_run_outputs(
    const std::vector< std::string > & whats
) const
{

    // Names of the viruses, given their ids
    auto vnames = [this](const std::vector< int > & ids) {
        std::vector< std::string > names;
        names.reserve(ids.size());
        for (auto i : ids)
            names.push_back(virus_name[i]);
        return names;
    };

    // Names of the layers, given their ids ("" for -1: not from a contact)
    auto lnames = [this](const std::vector< int > & ids) {
        std::vector< std::string > names;
        names.reserve(ids.size());
        for (auto i : ids)
            names.push_back(
                i < 0 ? std::string() :
                model->get_layer(static_cast< size_t >(i)).get_name()
            );
        return names;
    };

    // The status and stage columns of a vector of states
    auto add_pairs = [this](
        OutputTable & t,
        const std::vector< epiworld_fast_uint > & states,
        const std::string & prefix
    ) {
        std::vector< std::string > status, stage;
        split_states(states, &status, &stage);
        t.add(prefix + "status", status);
        t.add(prefix + "stage", stage);
    };

    RunOutputs out;

    for (const auto & what : whats)
    {

        OutputTable t;

        if (what == "virus_info")
        {
            std::vector< int > id(virus_name.size());
            std::iota(id.begin(), id.end(), 0);
            t.add("virus_id", id);
            t.add("virus", virus_name);
            t.add("date_recorded", virus_origin_date);
            t.add("parent", virus_parent_id);
        }
        else if (what == "virus_hist")
        {
            t.add("date", hist_virus_date);
            t.add("virus_id", hist_virus_id);
            t.add("virus", vnames(hist_virus_id));
            add_pairs(t, hist_virus_state, "");
            t.add("n", hist_virus_counts);
        }
        else if (what == "tool_info")
        {
            std::vector< int > id(tool_name.size());
            std::iota(id.begin(), id.end(), 0);
            t.add("id", id);
            t.add("tool_name", tool_name);
            t.add("date_recorded", tool_origin_date);
        }
        else if (what == "tool_hist")
        {
            t.add("date", hist_tool_date);
            t.add("id", hist_tool_id);
            add_pairs(t, hist_tool_state, "");
            t.add("n", hist_tool_counts);
        }
        else if (what == "total_hist")
        {
            t.add("date", hist_total_date);
            t.add("nviruses", hist_total_nviruses_active);
            add_pairs(t, hist_total_state, "");
            t.add("counts", hist_total_counts);
        }
        else if (what == "transmission")
        {
            t.add("date", transmission_date);
            t.add("virus_id", transmission_virus);
            t.add("virus", vnames(transmission_virus));
            t.add("source_exposure_date", transmission_source_exposure_date);
            t.add("source", transmission_source);
            t.add("target", transmission_target);
            t.add("layer_id", transmission_layer);
            t.add("layer", lnames(transmission_layer));
        }
        else if (what == "transition")
        {
            std::vector< int > date, counts;
            std::vector< epiworld_fast_uint > from, to;
            int ns = model->nstates;
            for (int i = 0; i <= model->today(); ++i)
                for (int f = 0; f < ns; ++f)
                    for (int s = 0; s < ns; ++s)
                    {
                        int n = hist_transition_matrix[i * (ns * ns) + s * ns + f];
                        if (n == 0)
                            continue;
                        date.push_back(i);
                        from.push_back(static_cast< epiworld_fast_uint >(f));
                        to.push_back(static_cast< epiworld_fast_uint >(s));
                        counts.push_back(n);
                    }
            t.add("date", date);
            add_pairs(t, from, "from_");
            add_pairs(t, to, "to_");
            t.add("counts", counts);
        }
        else if (what == "reproductive")
        {
            std::vector< int > id, source, date, rt;
            for (const auto & m : get_reproductive_number())
            {
                id.push_back(m.first[0u]);
                source.push_back(m.first[1u]);
                date.push_back(m.first[2u]);
                rt.push_back(m.second);
            }
            t.add("virus_id", id);
            t.add("virus", vnames(id));
            t.add("source", source);
            t.add("source_exposure_date", date);
            t.add("rt", rt);
        }
        else if (what == "generation")
        {
            std::vector< int > agent, virus, time, gentime;
            get_generation_time(agent, virus, time, gentime);
            t.add("virus", virus);
            t.add("source", agent);
            t.add("source_exposure_date", time);
            t.add("gentime", gentime);
        }
        else if (what == "active_cases" || what == "outbreak_size")
        {
            std::vector< int > date, id, n;
            if (what == "active_cases")
                get_active_cases(date, id, n);
            else
                get_outbreak_size(date, id, n);

            // Skipping the zeros
            std::vector< int > date_, id_, n_;
            for (size_t i = 0u; i < n.size(); ++i)
                if (n[i] > 0)
                {
                    date_.push_back(date[i]);
                    id_.push_back(id[i]);
                    n_.push_back(n[i]);
                }

            t.add("date", date_);
            t.add("virus_id", id_);
            t.add("virus", vnames(id_));
            t.add(what, n_);
        }
        else if (what == "hospitalizations")
        {
            std::vector< int > date, id, tool, count;
            std::vector< double > weight;
            get_hospitalizations(date, id, tool, count, weight);

            // Skipping the zeros
            std::vector< int > date_, id_, tool_, count_;
            std::vector< double > weight_;
            for (size_t i = 0u; i < date.size(); ++i)
                if (count[i] > 0 || weight[i] > 0.0)
                {
                    date_.push_back(date[i]);
                    id_.push_back(id[i]);
                    tool_.push_back(tool[i]);
                    count_.push_back(count[i]);
                    weight_.push_back(weight[i]);
                }

            t.add("date", date_);
            t.add("virus_id", id_);
            t.add("tool_id", tool_);
            t.add("count", count_);
            t.add("weight", weight_);
        }
        else
        {
            std::string valid;
            for (const auto & n : run_output_names())
                valid += " \"" + n + "\"";
            throw std::invalid_argument(
                "Unknown output \"" + what + "\". Valid outputs are:" + valid + "."
            );
        }

        out[what] = std::move(t);

    }

    return out;

}

inline void DataBase::write_data(
    std::string fn_virus_info,
    std::string fn_virus_hist,
    std::string fn_tool_info,
    std::string fn_tool_hist,
    std::string fn_total_hist,
    std::string fn_transmission,
    std::string fn_transition,
    std::string fn_reproductive_number,
    std::string fn_generation_time,
    std::string fn_active_cases,
    std::string fn_outbreak_size,
    std::string fn_hospitalizations
) const
{

    // Same order as run_output_names()
    const std::vector< std::string > fns = {
        fn_virus_info, fn_virus_hist, fn_tool_info, fn_tool_hist,
        fn_total_hist, fn_transmission, fn_transition,
        fn_reproductive_number, fn_generation_time, fn_active_cases,
        fn_outbreak_size, fn_hospitalizations
    };

    std::vector< std::string > whats;
    for (size_t i = 0u; i < fns.size(); ++i)
        if (fns[i] != "")
            whats.push_back(run_output_names()[i]);

    auto out = get_run_outputs(whats);

    for (size_t i = 0u; i < fns.size(); ++i)
        if (fns[i] != "")
            write_table(fns[i], out.at(run_output_names()[i]));

}

inline void DataBase::record_transmission(
    int i,
    int j,
    int virus,
    int i_expo_date,
    int layer
) {

    transmission_date.push_back(model->today());
    transmission_source.push_back(i);
    transmission_target.push_back(j);
    transmission_virus.push_back(virus);
    transmission_source_exposure_date.push_back(i_expo_date);
    transmission_layer.push_back(layer);

}

inline size_t DataBase::get_n_viruses() const
{
    return virus_name.size();
}

inline size_t DataBase::get_n_tools() const
{
    return tool_name.size();
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

inline MapVec_type<int,int> DataBase::get_reproductive_number()
const {

    // Checking size
    MapVec_type<int,int> map;

    // Number of digits of maxid
    for (size_t i = 0u; i < transmission_date.size(); ++i)
    {
        // Fabricating id
        std::vector< int > h = {
            transmission_virus[i],
            transmission_source[i],
            transmission_source_exposure_date[i]
        };

        // Adding to counter
        if (map.find(h) == map.end())
            map[h] = 1;
        else
            map[h]++;

        // The target is added
        std::vector< int > h_target = {
            transmission_virus[i],
            transmission_target[i],
            transmission_date[i]
        };

        map[h_target] = 0;
        
    }

    return map;

}

inline void DataBase::get_reproductive_number(
    std::string fn
) const {

    write_table(fn, get_run_outputs({"reproductive"}).at("reproductive"));

}

inline std::vector< epiworld_double > DataBase::get_transition_probability(
    bool print,
    bool normalize
) const {

    const auto& states_labels = model->states_labels;
    size_t n_state = states_labels.size();
    size_t n_days   = model->get_ndays();
    std::vector< epiworld_double > res(n_state * n_state, 0.0);
    std::vector< epiworld_double > rowsums(n_state, 0.0);

    for (size_t t = 0; t < n_days; ++t)
    {

        for (size_t s_i = 0; s_i < n_state; ++s_i)
        {

            for (size_t s_j = 0u; s_j < n_state; ++s_j)
            {
                res[s_i + s_j * n_state] += (
                    hist_transition_matrix[
                        s_i + s_j * n_state +
                        t * (n_state * n_state)
                    ]
                );
                
                rowsums[s_i] += hist_transition_matrix[
                    s_i + s_j * n_state +
                    t * (n_state * n_state)
                ];
            
            }

        }

    }

    if (normalize)
    {
        for (size_t s_i = 0; s_i < n_state; ++s_i)
        {
            // Nothing to normalize if the row is zero
            if (rowsums[s_i] == 0)
                continue;

            for (size_t s_j = 0; s_j < n_state; ++s_j)
                res[s_i + s_j * n_state] /= rowsums[s_i];
        }
    }

    if (print)
    {   

        size_t nchar = 0u;
        for (auto & l : states_labels)
            if (l.length() > nchar)
                nchar = l.length();

        std::string fmt = " - %-" + std::to_string(nchar) + "s";

        std::string fmt_entry = " % 4.2f";
        if (!normalize)
        {
            nchar = 0u;
            for (auto & l: res)
            {
                std::string tmp = std::to_string(l);
                if (tmp.length() > nchar)
                    nchar = tmp.length();
            }

            fmt_entry = " % " + std::to_string(nchar) + ".0f";
        } 

        
        printf_epiworld("\nTransition Probabilities:\n");
        for (size_t s_i = 0u; s_i < n_state; ++s_i)
        {
            printf_epiworld(fmt.c_str(), states_labels[s_i].c_str());
            for (size_t s_j = 0u; s_j < n_state; ++s_j)
            {
                if (
                    std::isnan(res[s_i + s_j * n_state]) ||
                    (res[s_i + s_j * n_state] < 1e-10)
                )
                {
                    printf_epiworld("     -");
                } else {
                    printf_epiworld(
                        fmt_entry.c_str(), res[s_i + s_j * n_state]
                    );
                }
            }
            printf_epiworld("\n");
        }

        printf_epiworld("\n");

    }

    return res;


} 

#define VECT_MATCH(a, b, c) \
    EPI_DEBUG_FAIL_AT_TRUE(a.size() != b.size(), c) \
    for (size_t __i = 0u; __i < a.size(); ++__i) \
    {\
        EPI_DEBUG_FAIL_AT_TRUE(a[__i] != b[__i], c) \
    }

inline bool DataBase::operator==(const DataBase & other) const
{
    VECT_MATCH(
        virus_name,
        other.virus_name,
        "DataBase:: virus_name[i] don't match"
    )

    VECT_MATCH(
        virus_origin_date,
        other.virus_origin_date,
        "DataBase:: virus_origin_date[i] don't match"
    )

    VECT_MATCH(
        virus_parent_id,
        other.virus_parent_id,
        "DataBase:: virus_parent_id[i] don't match"
    )

    VECT_MATCH(
        tool_name,
        other.tool_name,
        "DataBase:: tool_name[i] don't match"
    )

    VECT_MATCH(
        tool_origin_date,
        other.tool_origin_date,
        "DataBase:: tool_origin_date[i] don't match"
    )

    
    EPI_DEBUG_FAIL_AT_TRUE(
        sampling_freq != other.sampling_freq,
        "DataBase:: sampling_freq don't match."
    )

    // Variants history
    VECT_MATCH(
        hist_virus_date,
        other.hist_virus_date,
        "DataBase:: hist_virus_date[i] don't match"
    )

    VECT_MATCH(
        hist_virus_id,
        other.hist_virus_id,
        "DataBase:: hist_virus_id[i] don't match"
    )

    VECT_MATCH(
        hist_virus_state,
        other.hist_virus_state,
        "DataBase:: hist_virus_state[i] don't match"
    )

    VECT_MATCH(
        hist_virus_counts,
        other.hist_virus_counts,
        "DataBase:: hist_virus_counts[i] don't match"
    )

    // Tools history
    VECT_MATCH(
        hist_tool_date,
        other.hist_tool_date,
        "DataBase:: hist_tool_date[i] don't match"
    )

    VECT_MATCH(
        hist_tool_id,
        other.hist_tool_id,
        "DataBase:: hist_tool_id[i] don't match"
    )

    VECT_MATCH(
        hist_tool_state,
        other.hist_tool_state,
        "DataBase:: hist_tool_state[i] don't match"
    )

    VECT_MATCH(
        hist_tool_counts,
        other.hist_tool_counts,
        "DataBase:: hist_tool_counts[i] don't match"
    )

    // Overall hist
    VECT_MATCH(
        hist_total_date,
        other.hist_total_date,
        "DataBase:: hist_total_date[i] don't match"
    )

    VECT_MATCH(
        hist_total_nviruses_active,
        other.hist_total_nviruses_active,
        "DataBase:: hist_total_nviruses_active[i] don't match"
    )

    VECT_MATCH(
        hist_total_state,
        other.hist_total_state,
        "DataBase:: hist_total_state[i] don't match"
    )

    VECT_MATCH(
        hist_total_counts,
        other.hist_total_counts,
        "DataBase:: hist_total_counts[i] don't match"
    )

    VECT_MATCH(
        hist_transition_matrix,
        other.hist_transition_matrix,
        "DataBase:: hist_transition_matrix[i] don't match"
    )

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    EPI_DEBUG_FAIL_AT_TRUE(
        today_virus.size() != other.today_virus.size(),
        "DataBase:: today_virus.size() don't match."
    )
    
    for (size_t i = 0u; i < today_virus.size(); ++i)
    {
        VECT_MATCH(
            today_virus[i], other.today_virus[i],
            "DataBase:: today_virus[i] don't match"
            )
    }

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    EPI_DEBUG_FAIL_AT_TRUE(
        today_tool.size() != other.today_tool.size(),
        "DataBase:: today_tool.size() don't match."
    )
    
    for (size_t i = 0u; i < today_tool.size(); ++i)
    {
        VECT_MATCH(
            today_tool[i], other.today_tool[i],
            "DataBase:: today_tool[i] don't match"
            )
    }

    // {Susceptible, Infected, etc.}
    VECT_MATCH(
        today_total, other.today_total,
        "DataBase:: today_total[i] don't match"
        )

    // Totals
    EPI_DEBUG_FAIL_AT_TRUE(
        today_total_nviruses_active != other.today_total_nviruses_active,
        "DataBase:: today_total_nviruses_active don't match."
    )

    // Transmission network
    VECT_MATCH( ///< Date of the transmission eve
        transmission_date,
        other.transmission_date,
        "DataBase:: transmission_date[i] don't match"
    )

    VECT_MATCH( ///< Id of the sour
        transmission_source,
        other.transmission_source,
        "DataBase:: transmission_source[i] don't match"
    )

    VECT_MATCH( ///< Id of the targ
        transmission_target,
        other.transmission_target,
        "DataBase:: transmission_target[i] don't match"
    )

    VECT_MATCH( ///< Id of the varia
        transmission_virus,
        other.transmission_virus,
        "DataBase:: transmission_virus[i] don't match"
    )

    VECT_MATCH( ///< Date when the source acquired the varia
        transmission_source_exposure_date,
        other.transmission_source_exposure_date,
        "DataBase:: transmission_source_exposure_date[i] don't match"
    )

    VECT_MATCH(
        transmission_layer,
        other.transmission_layer,
        "DataBase:: transmission_layer[i] don't match"
    )

    VECT_MATCH(
        transition_matrix,
        other.transition_matrix,
        "DataBase:: transition_matrix[i] don't match"
    )

    return true;

}

inline void DataBase::get_generation_time(
    std::vector< int > & agent_id,
    std::vector< int > & virus_id,
    std::vector< int > & time,
    std::vector< int > & gentime
) const {
    
    size_t nevents = transmission_date.size();

    agent_id.reserve(nevents);
    virus_id.reserve(nevents);
    time.reserve(nevents);
    gentime.reserve(nevents);

    // Iterating through the individuals
    for (size_t i = 0u; i < nevents; ++i)
    {
        int agent_id_i = transmission_target[i];
        agent_id.push_back(agent_id_i);
        virus_id.push_back(transmission_virus[i]);
        time.push_back(transmission_date[i]);

        bool found = false;
        for (size_t j = i; j < nevents; ++j)
        {

            if (transmission_source[j] == agent_id_i)
            {
                gentime.push_back(transmission_date[j] - time[i]);
                found = true;
                break;
            }

        }

        // If there's no transmission, we set the generation time to
        // minus 1;
        if (!found)
            gentime.push_back(-1);

    }

    agent_id.shrink_to_fit();
    virus_id.shrink_to_fit();
    time.shrink_to_fit();
    gentime.shrink_to_fit();

    return;

}

inline void DataBase::get_generation_time(
    std::string fn
) const
{

    write_table(fn, get_run_outputs({"generation"}).at("generation"));

}

inline void DataBase::record_hospitalization(Agent & agent)
{
    m_hospitalizations.record(agent, *model);
}

inline void DataBase::get_hospitalizations(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & tool_id,
    std::vector<int> & count,
    std::vector<double> & weight
) const
{
    // Get the number of days from the model and add 1 since days are 0-indexed
    // today() returns the last day, so we need today() + 1 for the full count
    int ndays = model->today() + 1;
    m_hospitalizations.get(ndays, date, virus_id, tool_id, count, weight);
}

#undef VECT_MATCH

inline void DataBase::get_hist_transitions(
    std::vector< int > & date,
    std::vector< std::string > & transition,
    std::vector< int > & counts
) const
{

    if (!model->nh)
        throw std::logic_error(
            "get_hist_transitions() needs a model with a natural history."
        );

    date = hist_transition_date;
    counts = hist_transition_counts;
    // Stage transitions first, then the status ones
    const int n_stage = static_cast< int >(model->nh->get_n_transitions());
    transition.resize(hist_transition_id.size());
    for (size_t k = 0u; k < hist_transition_id.size(); ++k)
    {
        const int t = hist_transition_id[k];
        transition[k] = (t < n_stage) ?
            model->nh->get_transition(t).name :
            model->st->get_transition(t - n_stage).name;
    }

}

inline void SaverMemory::operator()(size_t sim_id, Model * model) const
{
    auto out = model->get_db().get_run_outputs(state->whats);
    if (sim_id >= state->runs.size())
        state->runs.resize(sim_id + 1u);
    state->runs[sim_id] = std::move(out);
}

#endif
