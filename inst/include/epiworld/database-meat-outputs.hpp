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

    // The status and stage columns of a vector of cells of disease `d`
    auto add_pairs = [this](
        OutputTable & t,
        int d,
        const std::vector< epiworld_fast_uint > & cells,
        const std::string & prefix
    ) {
        std::vector< std::string > status, stage;
        split_states(d, cells, &status, &stage);
        t.add(prefix + "status", status);
        t.add(prefix + "stage", stage);
    };

    // Keeps the rows of disease `d`'s viruses (with several diseases, the
    // per-virus outputs list every virus in every disease's computation)
    auto keep_virus_rows = [this](
        int d,
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< int > & n
    ) {
        std::vector< int > date_, id_, n_;
        for (size_t i = 0u; i < n.size(); ++i)
            if ((n[i] > 0) && ((recs.size() == 1u) || (disease_of_virus(id[i]) == d)))
            {
                date_.push_back(date[i]);
                id_.push_back(id[i]);
                n_.push_back(n[i]);
            }
        date.swap(date_);
        id.swap(id_);
        n.swap(n_);
    };

    RunOutputs out;

    for (const auto & what : whats)
    {

        // The tables of the registry: once
        if (what == "virus_info")
        {
            OutputTable t;
            std::vector< int > id(reg.virus_name.size());
            std::iota(id.begin(), id.end(), 0);
            std::vector< std::string > names;
            for (int v : id)
                disease_rows(&names, 1u, model->tracks.empty() ? 0 : disease_of_virus(v));
            t.add("disease", names);
            t.add("virus_id", id);
            t.add("virus", reg.virus_name);
            t.add("date_recorded", reg.virus_origin_date);
            t.add("parent", reg.virus_parent_id);
            out[what] = std::move(t);
            continue;
        }

        if (what == "tool_info")
        {
            OutputTable t;
            std::vector< int > id(reg.tool_name.size());
            std::iota(id.begin(), id.end(), 0);
            t.add("id", id);
            t.add("tool_name", reg.tool_name);
            t.add("date_recorded", reg.tool_origin_date);
            out[what] = std::move(t);
            continue;
        }

        if (
            (what != "virus_hist") && (what != "tool_hist") &&
            (what != "total_hist") && (what != "transmission") &&
            (what != "transition") && (what != "reproductive") &&
            (what != "generation") && (what != "active_cases") &&
            (what != "outbreak_size") && (what != "hospitalizations")
        )
        {
            std::string valid;
            for (const auto & n : run_output_names())
                valid += " \"" + n + "\"";
            throw std::invalid_argument(
                "Unknown output \"" + what + "\". Valid outputs are:" + valid + "."
            );
        }

        // The disease's tables, one disease after the other, with the disease
        // of each row first
        OutputTable all;
        for (size_t di = 0u; di < recs.size(); ++di)
        {

            const int d = static_cast< int >(di);
            const auto & r = recs[di];
            OutputTable t;

            if (what == "virus_hist")
            {
                t.add("date", r.hist_virus_date);
                t.add("virus_id", r.hist_virus_id);
                t.add("virus", vnames(r.hist_virus_id));
                add_pairs(t, d, r.hist_virus_state, "");
                t.add("n", r.hist_virus_counts);
            }
            else if (what == "tool_hist")
            {
                t.add("date", r.hist_tool_date);
                t.add("id", r.hist_tool_id);
                add_pairs(t, d, r.hist_tool_state, "");
                t.add("n", r.hist_tool_counts);
            }
            else if (what == "total_hist")
            {
                t.add("date", r.hist_total_date);
                t.add("nviruses", r.hist_total_nviruses_active);
                add_pairs(t, d, r.hist_total_state, "");
                t.add("counts", r.hist_total_counts);
            }
            else if (what == "transmission")
            {
                t.add("date", r.transmission_date);
                t.add("virus_id", r.transmission_virus);
                t.add("virus", vnames(r.transmission_virus));
                t.add("source_exposure_date", r.transmission_source_exposure_date);
                t.add("source", r.transmission_source);
                t.add("target", r.transmission_target);
                t.add("layer_id", r.transmission_layer);
                t.add("layer", lnames(r.transmission_layer));
            }
            else if (what == "transition")
            {
                std::vector< int > date, counts;
                std::vector< epiworld_fast_uint > from, to;
                const int ns = static_cast< int >(r.ncells);
                for (int i = 0; i <= model->today(); ++i)
                    for (int f = 0; f < ns; ++f)
                        for (int s = 0; s < ns; ++s)
                        {
                            int n = r.hist_transition_matrix[i * (ns * ns) + s * ns + f];
                            if (n == 0)
                                continue;
                            date.push_back(i);
                            from.push_back(static_cast< epiworld_fast_uint >(f));
                            to.push_back(static_cast< epiworld_fast_uint >(s));
                            counts.push_back(n);
                        }
                t.add("date", date);
                add_pairs(t, d, from, "from_");
                add_pairs(t, d, to, "to_");
                t.add("counts", counts);
            }
            else if (what == "reproductive")
            {
                std::vector< int > id, source, date, rt;
                for (const auto & m : get_reproductive_number(model->tracks.empty() ? all_diseases : d))
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
                get_generation_time(agent, virus, time, gentime, model->tracks.empty() ? all_diseases : d);
                t.add("virus", virus);
                t.add("source", agent);
                t.add("source_exposure_date", time);
                t.add("gentime", gentime);
            }
            else if (what == "active_cases" || what == "outbreak_size")
            {
                std::vector< int > date, id, n;
                if (what == "active_cases")
                    get_active_cases(date, id, n, model->tracks.empty() ? all_diseases : d);
                else
                    get_outbreak_size(date, id, n, model->tracks.empty() ? all_diseases : d);

                // Skipping the zeros
                keep_virus_rows(d, date, id, n);

                t.add("date", date);
                t.add("virus_id", id);
                t.add("virus", vnames(id));
                t.add(what, n);
            }
            else if (what == "hospitalizations")
            {
                std::vector< int > date, id, tool, count;
                std::vector< double > weight;
                get_hospitalizations(date, id, tool, count, weight, model->tracks.empty() ? all_diseases : d);

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

            std::vector< std::string > names;
            disease_rows(&names, t.nrow(), d);
            t.colnames.insert(t.colnames.begin(), "disease");
            t.columns.insert(t.columns.begin(), std::move(names));

            all.append(t);

        }

        out[what] = std::move(all);

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
