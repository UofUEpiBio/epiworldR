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

inline size_t DataBase::ncells_of(int d) const
{
    if (model->tracks.empty())
        return static_cast< size_t >(model->nstates);
    return model->nh_nstatus() * model->tracks[static_cast< size_t >(d)].n;
}

inline epiworld_fast_uint DataBase::cell_of(int d, epiworld_fast_uint state) const
{
    // With one disease (or none), the engine's state is the cell
    if (model->tracks.size() <= 1u)
        return state;
    const unsigned int c = static_cast< unsigned int >(state);
    return static_cast< epiworld_fast_uint >(
        static_cast< unsigned int >(model->nh_status_of(c)) +
        static_cast< unsigned int >(model->nh_nstatus()) * static_cast< unsigned int >(
            model->nh_stage_of(model->tracks[static_cast< size_t >(d)], c)
        )
    );
}

inline std::vector< int > DataBase::diseases_of(int disease) const
{
    check_disease(disease);
    if (disease != all_diseases)
        return {disease};
    std::vector< int > ds(recs.size());
    std::iota(ds.begin(), ds.end(), 0);
    return ds;
}

inline int DataBase::disease_of_virus(int virus_id) const
{
    const int d = reg.virus_disease[static_cast< size_t >(virus_id)];
    return d < 0 ? 0 : d;
}

inline void DataBase::count_virus(int virus_id, epiworld_fast_uint state, int delta)
{
    const int d = (recs.size() == 1u) ? 0 : disease_of_virus(virus_id);
    auto & r = recs[static_cast< size_t >(d)];
    #ifdef EPI_DEBUG
    r.today_virus.at(static_cast< size_t >(virus_id)).at(cell_of(d, state)) += delta;
    #else
    r.today_virus[static_cast< size_t >(virus_id)][cell_of(d, state)] += delta;
    #endif
}

inline void DataBase::count_tool(int tool_id, epiworld_fast_uint state, int delta)
{
    // One disease: a tool always acts on it (see DiseaseTargets::resolve())
    if (recs.size() == 1u)
    {
        recs[0u].today_tool[static_cast< size_t >(tool_id)][state] += delta;
        return;
    }

    const uint64_t mask = reg.tool_mask[static_cast< size_t >(tool_id)];
    for (size_t d = 0u; d < recs.size(); ++d)
        if ((mask >> d) & 1u)
            recs[d].today_tool[static_cast< size_t >(tool_id)][
                cell_of(static_cast< int >(d), state)
            ] += delta;
}

inline void DataBase::reset()
{

    // One record per disease (one if the model has none)
    recs.resize(std::max< size_t >(1u, model->tracks.size()));

    // The last run's variants go; the founders keep their ids
    drop_variants();

    for (size_t di = 0u; di < recs.size(); ++di)
    {

        const int d = static_cast< int >(di);
        auto & r = recs[di];
        const size_t nc = ncells_of(d);
        r.ncells = nc;

        // Initializing the counts
        r.today_total.assign(nc, 0);
        for (auto & p : model->get_agents())
            ++r.today_total[cell_of(d, p.get_state())];

        #ifdef EPI_DEBUG
        // Only the first should be different from zero
        {
            int n = static_cast<int>(model->size());
            if (r.today_total[0] != n)
                throw std::runtime_error("The number of susceptible agents is not equal to the total number of agents.");

            if (std::accumulate(r.today_total.begin(), r.today_total.end(), 0) != n)
                throw std::runtime_error("The total number of agents is not equal to the sum of the number of agents in each state.");
        }
        #endif

        r.transition_matrix.assign(nc * nc, 0);
        for (size_t s = 0u; s < nc; ++s)
            r.transition_matrix[s + s * nc] = r.today_total[s];

        r.hist_virus_date.clear();
        r.hist_virus_id.clear();
        r.hist_virus_state.clear();
        r.hist_virus_counts.clear();

        r.hist_tool_date.clear();
        r.hist_tool_id.clear();
        r.hist_tool_state.clear();
        r.hist_tool_counts.clear();

        r.today_virus.resize(get_n_viruses());
        for (auto & virus_states : r.today_virus)
            virus_states.assign(nc, 0);

        r.today_tool.resize(get_n_tools());
        for (auto & tool_states : r.today_tool)
            tool_states.assign(nc, 0);

        // The disease's variants
        r.today_total_nviruses_active = 0;
        for (size_t v = 0u; v < reg.virus_disease.size(); ++v)
            if (disease_of_virus(static_cast< int >(v)) == d)
                ++r.today_total_nviruses_active;

        r.hist_total_date.clear();
        r.hist_total_state.clear();
        r.hist_total_nviruses_active.clear();
        r.hist_total_counts.clear();
        r.hist_transition_matrix.clear();

        r.transmission_date.clear();
        r.transmission_virus.clear();
        r.transmission_source.clear();
        r.transmission_target.clear();
        r.transmission_source_exposure_date.clear();
        r.transmission_layer.clear();

        // The disease's declared transitions (the statuses', once: after the
        // first disease's)
        r.transition_today.assign(
            !model->tracks.empty() ?
                model->tracks[di].nh->get_n_transitions() +
                    (((di == 0u) && model->st) ? model->st->get_n_transitions() : 0u) :
                0u,
            0
        );
        r.hist_transition_date.clear();
        r.hist_transition_id.clear();
        r.hist_transition_counts.clear();

        // By variant: the transitions an agent with the virus can take
        r.virus_transitions.clear();
        if (!model->tracks.empty())
        {
            const auto & nh = *model->tracks[di].nh;
            for (size_t t = 0u; t < nh.get_n_transitions(); ++t)
            {
                const auto & tr = nh.get_transition(static_cast< int >(t));
                bool carries = (tr.kind == TransitionKind::infection);
                for (int s : tr.from)
                    carries = carries || nh.is_carrier(s);
                if (carries)
                    r.virus_transitions.push_back(static_cast< int >(t));
            }
        }
        r.transition_virus_today.resize(get_n_viruses());
        for (auto & c : r.transition_virus_today)
            c.clear();
        r.hist_transition_virus_date.clear();
        r.hist_transition_virus_id.clear();
        r.hist_transition_virus_transition.clear();
        r.hist_transition_virus_counts.clear();

        r.m_hospitalizations.reset();

    }

    return;

}

inline DataBase::DataBase(const DataBase & db) :
    reg(db.reg),
    recs(db.recs),
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
    #ifdef EPI_DEBUG
    // Checking whether the sums correspond, in each disease
    for (size_t di = 0u; di < recs.size(); ++di)
    {
        const auto & r = recs[di];
        EPI_DEBUG_SUM_INT(r.today_total, model->size())
        EPI_DEBUG_ALL_NON_NEGATIVE(r.today_total)

        std::vector< int > _today_total_cp(r.today_total.size(), 0);
        for (auto & p : model->population)
            _today_total_cp[cell_of(static_cast< int >(di), p.get_state())]++;

        EPI_DEBUG_VECTOR_MATCH_INT(
            _today_total_cp, r.today_total,
            "Sums of __today_total_cp in database-meat.hpp"
            )
    }

    // EPI_DEBUG only: the model's agents-by-state index (see
    // Model::get_agents_in_state()) is maintained incrementally, event by
    // event, so check it against the population once per step. Each state's
    // block must hold exactly the agents in that state, each agent's recorded
    // position must be right, and the degree and carrier sums used by the
    // transmission step must match a recount. A drift here would silently
    // bias which agents push, or which mode "auto" picks.
    if (model->state_index_ready)
    {

        std::vector< size_t > _deg(model->nstates, 0u);
        std::vector< size_t > _carriers(model->nstates, 0u);
        std::vector< size_t > _carrier_deg(model->nstates, 0u);
        std::vector< size_t > _members(model->nstates, 0u);
        for (auto & p : model->population)
        {
            const size_t _d = model->network_degree(static_cast< size_t >(p.get_id()));
            _deg[p.get_state()] += _d;
            _members[p.get_state()]++;
            if (p.get_virus() != nullptr)
            {
                _carriers[p.get_state()]++;
                _carrier_deg[p.get_state()] += _d;
            }
        }

        for (size_t s = 0u; s < model->nstates; ++s)
        {

            const auto members = model->state_index_members(s);
            if (members.size() != _members[s])
                throw std::logic_error("[epi-debug] DataBase::record state index size doesn't match the population.");

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
        for (const auto & r : recs)
            if (
                !r.hist_total_date.empty() || !r.hist_total_nviruses_active.empty() ||
                !r.hist_total_state.empty() || !r.hist_total_counts.empty() ||
                !r.hist_virus_date.empty() || !r.hist_virus_id.empty() ||
                !r.hist_virus_state.empty() || !r.hist_virus_counts.empty() ||
                !r.hist_tool_date.empty() || !r.hist_tool_id.empty() ||
                !r.hist_tool_state.empty() || !r.hist_tool_counts.empty()
            )
                EPI_DEBUG_ERROR(std::logic_error, "DataBase::record the histories should be empty on day 0.")
    #endif
    ////////////////////////////////////////////////////////////////////////////

    // Only store every now and then
    if ((model->today() % sampling_freq) == 0)
        for (auto & r : recs)
            record_disease(r);

}

inline void DataBase::record_disease(DiseaseRecord & r)
{

    const size_t nc = r.ncells;

    // Recording virus's history (the disease's variants)
    const int d = static_cast< int >(&r - recs.data());
    for (size_t v = 0u; v < r.today_virus.size(); ++v)
    {

        if ((recs.size() > 1u) && (disease_of_virus(static_cast< int >(v)) != d))
            continue;

        for (epiworld_fast_uint s = 0u; s < nc; ++s)
        {

            r.hist_virus_date.push_back(model->today());
            r.hist_virus_id.push_back(static_cast< int >(v));
            r.hist_virus_state.push_back(s);

            #ifdef EPI_DEBUG
            if (r.today_virus[v][s] < 0)
                throw std::logic_error(
                    "The count of viruses in DataBase::record cannot be negative."
                );
            #endif

            r.hist_virus_counts.push_back(r.today_virus[v][s]);

        }

    }

    // Recording tool's history (in the diseases each targets)
    for (size_t t = 0u; t < r.today_tool.size(); ++t)
    {

        if (!((reg.tool_mask[t] >> d) & 1u))
            continue;

        for (epiworld_fast_uint s = 0u; s < nc; ++s)
        {

            r.hist_tool_date.push_back(model->today());
            r.hist_tool_id.push_back(static_cast< int >(t));
            r.hist_tool_state.push_back(s);
            r.hist_tool_counts.push_back(r.today_tool[t][s]);

        }

    }

    // Recording the overall history
    for (epiworld_fast_uint s = 0u; s < nc; ++s)
    {
        r.hist_total_date.push_back(model->today());
        r.hist_total_nviruses_active.push_back(r.today_total_nviruses_active);
        r.hist_total_state.push_back(s);
        r.hist_total_counts.push_back(r.today_total[s]);
    }

    for (const auto& cell : r.transition_matrix)
        r.hist_transition_matrix.push_back(cell);

    for (size_t k = 0u; k < r.transition_today.size(); ++k)
    {
        r.hist_transition_date.push_back(model->today());
        r.hist_transition_id.push_back(static_cast< int >(k));
        r.hist_transition_counts.push_back(r.transition_today[k]);
        r.transition_today[k] = 0;
    }

    // By variant (the disease's variants)
    if (!r.virus_transitions.empty())
        for (size_t v = 0u; v < r.transition_virus_today.size(); ++v)
        {

            if ((recs.size() > 1u) && (disease_of_virus(static_cast< int >(v)) != d))
                continue;

            auto & c = r.transition_virus_today[v];
            for (int t : r.virus_transitions)
            {
                r.hist_transition_virus_date.push_back(model->today());
                r.hist_transition_virus_id.push_back(static_cast< int >(v));
                r.hist_transition_virus_transition.push_back(t);
                r.hist_transition_virus_counts.push_back(
                    c.empty() ? 0 : c[static_cast< size_t >(t)]
                );
            }
            std::fill(c.begin(), c.end(), 0);

        }

    // Now the diagonal must reflect the state
    for (size_t s_i = 0u; s_i < nc; ++s_i)
    {

        for (size_t s_j = 0u; s_j < nc; ++s_j)
        {

            if ((s_i != s_j) && (r.transition_matrix[s_i + s_j * nc] > 0))
            {
                r.transition_matrix[s_j + s_j * nc] +=
                    r.transition_matrix[s_i + s_j * nc];

                r.transition_matrix[s_i + s_j * nc] = 0;
            }

        }

    }

    #ifdef EPI_DEBUG
    for (size_t s_i = 0u; s_i < nc; ++s_i)
    {
        if (r.transition_matrix[s_i + s_i * nc] != r.today_total[s_i])
            throw std::logic_error(
                "The diagonal of the updated transition Matrix should match the daily totals"
                );
    }
    #endif

}

inline void DataBase::drop_variants()
{

    const size_t n = reg.n_founders;
    if (reg.virus_name.size() <= n)
        return;

    reg.virus_name.resize(n);
    reg.virus_origin_date.resize(n);
    reg.virus_parent_id.resize(n);
    reg.virus_disease.resize(n);
    for (auto & r : recs)
    {
        r.today_virus.resize(n);
        r.transition_virus_today.resize(n);
    }

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
    reg.virus_disease.push_back(v.disease);

    if (recs.empty())
        recs.resize(1u);

    for (size_t d = 0u; d < recs.size(); ++d)
    {
        recs[d].today_virus.push_back(
            std::vector< int >(recs[d].ncells, 0)
        );
        recs[d].transition_virus_today.emplace_back();
    }
    if (static_cast< size_t >(disease_of_virus(id)) < recs.size())
        recs[static_cast< size_t >(disease_of_virus(id))].today_total_nviruses_active++;

    v.set_id(id);
    v.set_date(model->today());

    // A founder starts a lineage, which its mutations keep
    if (founder)
        v.lineage_id = id;
    else if (v.get_agent() != nullptr)
    {
        // The agent's counts move from the parent to the variant
        const auto s = v.get_agent()->get_state();
        count_virus(parent, s, -1);
        count_virus(id, s, 1);
    }

}

inline void DataBase::record_tool(Tool & t)
{

    t.set_id(static_cast< int >(reg.tool_name.size()));
    t.set_date(model->today());
    reg.tool_name.push_back(t.get_name());
    reg.tool_origin_date.push_back(model->today());
    reg.tool_mask.push_back(t.targets_.mask);

    if (recs.empty())
        recs.resize(1u);

    for (auto & r : recs)
        r.today_tool.push_back(std::vector< int >(r.ncells, 0));

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
    #endif

    // One disease: the engine's states are its cells
    if (recs.size() == 1u)
    {
        auto & r = recs[0u];
        if (undo)
        {
            r.today_total[prev_state]++;
            r.today_total[new_state]--;
        }
        else
        {
            r.today_total[prev_state]--;
            r.today_total[new_state]++;
        }
        record_transition(r, prev_state, new_state, undo);
        return;
    }

    update_state_several(prev_state, new_state, undo);

}

inline void DataBase::update_state_several(
    epiworld_fast_uint prev_state,
    epiworld_fast_uint new_state,
    bool undo
)
{

    // Each disease counts its own cells (a status change moves the agent
    // in every disease; a stage change, in one)
    for (size_t di = 0u; di < recs.size(); ++di)
    {

        const int d = static_cast< int >(di);
        const epiworld_fast_uint from = cell_of(d, prev_state);
        const epiworld_fast_uint to = cell_of(d, new_state);
        if (from == to)
            continue;

        auto & r = recs[di];
        if (undo)
        {
            r.today_total[from]++;
            r.today_total[to]--;
        }
        else
        {
            r.today_total[from]--;
            r.today_total[to]++;
        }

        record_transition(r, from, to, undo);

    }

}

inline void DataBase::update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    const int sign = undo ? -1 : 1;
    count_virus(static_cast< int >(virus_id), prev_state, -sign);
    count_virus(static_cast< int >(virus_id), new_state, sign);

    #ifdef EPI_DEBUG
    const int d = disease_of_virus(static_cast< int >(virus_id));
    const auto & tv = recs[static_cast< size_t >(d)].today_virus[virus_id];
    if ((tv[cell_of(d, prev_state)] < 0) || (tv[cell_of(d, new_state)] < 0))
        throw std::logic_error("An entry in today_virus is negative.");
    #endif

    return;
    
}

inline void DataBase::update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

    const int sign = undo ? -1 : 1;
    count_tool(static_cast< int >(tool_id), prev_state, -sign);
    count_tool(static_cast< int >(tool_id), new_state, sign);

    return;

}

inline void DataBase::record_transition(
    DiseaseRecord & r,
    epiworld_fast_uint from,
    epiworld_fast_uint to,
    bool undo
) {

    const size_t nc = r.ncells;
    if (undo)
    {

        r.transition_matrix[to * nc + from]--;
        r.transition_matrix[from * nc + from]++;

    } else {

        r.transition_matrix[to * nc + from]++;
        r.transition_matrix[from * nc + from]--;

    }

    #ifdef EPI_DEBUG
    if (r.transition_matrix[from * nc + from] < 0)
        throw std::logic_error("An entry in transition matrix is negative.");
    #endif

}

inline void DataBase::record_declared_transition(int track, int t, int virus)
{

    if (Model::nh_is_status_tr(t))
    {
        ++recs[0u].transition_today[
            model->tracks[0u].nh->get_n_transitions() +
            static_cast< size_t >(Model::nh_status_tr(t))
        ];
        return;
    }

    auto & r = recs[static_cast< size_t >(track)];
    ++r.transition_today[static_cast< size_t >(t)];

    // By variant, if the agent carries (or caught) one of the disease
    if ((virus < 0) || ((recs.size() > 1u) && (disease_of_virus(virus) != track)))
        return;

    auto & c = r.transition_virus_today[static_cast< size_t >(virus)];
    if (c.empty())
        c.assign(r.transition_today.size(), 0);
    ++c[static_cast< size_t >(t)];

}

inline void DataBase::record_transmission(
    int i,
    int j,
    int virus,
    int i_expo_date,
    int layer
) {

    auto & r = recs[static_cast< size_t >(disease_of_virus(virus))];
    r.transmission_date.push_back(model->today());
    r.transmission_source.push_back(i);
    r.transmission_target.push_back(j);
    r.transmission_virus.push_back(virus);
    r.transmission_source_exposure_date.push_back(i_expo_date);
    r.transmission_layer.push_back(layer);

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
    // Under the disease of the virus the agent carries
    const auto & v = agent.get_virus();
    const int d = (v != nullptr) ? disease_of_virus(v->get_id()) : 0;
    recs[static_cast< size_t >(d)].m_hospitalizations.record(agent, *model);
}

#endif
