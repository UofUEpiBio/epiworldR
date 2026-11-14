#ifndef EPIWORLD_DATABASE_MEAT_OUTPUTS_HPP
#define EPIWORLD_DATABASE_MEAT_OUTPUTS_HPP

/**
 * @file database-meat-outputs.hpp
 * @brief The run's output tables (`get_run_outputs()`, used by the savers)
 * and the files `write_data()` writes.
 */

inline RunOutputs DataBase::get_run_outputs(
    const std::vector< std::string > & whats
) const
{

    // Names of the viruses, given their ids
    auto vnames = [this](const std::vector< int > & ids) {
        std::vector< std::string > names;
        names.reserve(ids.size());
        for (auto i : ids)
            names.push_back(reg.virus_name[i]);
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
            std::vector< int > id(reg.virus_name.size());
            std::iota(id.begin(), id.end(), 0);
            t.add("virus_id", id);
            t.add("virus", reg.virus_name);
            t.add("date_recorded", reg.virus_origin_date);
            t.add("parent", reg.virus_parent_id);
        }
        else if (what == "virus_hist")
        {
            t.add("date", rec.hist_virus_date);
            t.add("virus_id", rec.hist_virus_id);
            t.add("virus", vnames(rec.hist_virus_id));
            add_pairs(t, rec.hist_virus_state, "");
            t.add("n", rec.hist_virus_counts);
        }
        else if (what == "tool_info")
        {
            std::vector< int > id(reg.tool_name.size());
            std::iota(id.begin(), id.end(), 0);
            t.add("id", id);
            t.add("tool_name", reg.tool_name);
            t.add("date_recorded", reg.tool_origin_date);
        }
        else if (what == "tool_hist")
        {
            t.add("date", rec.hist_tool_date);
            t.add("id", rec.hist_tool_id);
            add_pairs(t, rec.hist_tool_state, "");
            t.add("n", rec.hist_tool_counts);
        }
        else if (what == "total_hist")
        {
            t.add("date", rec.hist_total_date);
            t.add("nviruses", rec.hist_total_nviruses_active);
            add_pairs(t, rec.hist_total_state, "");
            t.add("counts", rec.hist_total_counts);
        }
        else if (what == "transmission")
        {
            t.add("date", rec.transmission_date);
            t.add("virus_id", rec.transmission_virus);
            t.add("virus", vnames(rec.transmission_virus));
            t.add("source_exposure_date", rec.transmission_source_exposure_date);
            t.add("source", rec.transmission_source);
            t.add("target", rec.transmission_target);
            t.add("layer_id", rec.transmission_layer);
            t.add("layer", lnames(rec.transmission_layer));
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
                        int n = rec.hist_transition_matrix[i * (ns * ns) + s * ns + f];
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

        // The disease of each row, first (tools are not one disease's)
        if (what != "tool_info")
        {
            std::vector< std::string > names;
            disease_rows(&names, t.nrow());
            t.colnames.insert(t.colnames.begin(), "disease");
            t.columns.insert(t.columns.begin(), std::move(names));
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

inline void DataBase::get_reproductive_number(
    std::string fn
) const {

    write_table(fn, get_run_outputs({"reproductive"}).at("reproductive"));

}

inline void DataBase::get_generation_time(
    std::string fn
) const
{

    write_table(fn, get_run_outputs({"generation"}).at("generation"));

}

inline void SaverMemory::operator()(size_t sim_id, Model * model) const
{
    auto out = model->get_db().get_run_outputs(state->whats);
    if (sim_id >= state->runs.size())
        state->runs.resize(sim_id + 1u);
    state->runs[sim_id] = std::move(out);
}

#endif
