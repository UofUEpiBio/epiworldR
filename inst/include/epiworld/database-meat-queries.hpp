#ifndef EPIWORLD_DATABASE_MEAT_QUERIES_HPP
#define EPIWORLD_DATABASE_MEAT_QUERIES_HPP

/**
 * @file database-meat-queries.hpp
 * @brief Reading what the database recorded: counts, histories,
 * transitions, transmissions, reproductive numbers and hospitalizations.
 * @details Each function that returns rows returns those of the diseases its
 * `disease` argument asks for (every disease by default), one disease after
 * the other, and their disease's name in `disease_names`.
 */

/// Appends `src` to `dst`
template< typename T >
inline void epi_append(std::vector< T > & dst, const std::vector< T > & src)
{
    dst.insert(dst.end(), src.begin(), src.end());
}

inline void DataBase::split_state(
    int d,
    epiworld_fast_uint cell,
    std::string & status,
    std::string & stage
) const
{

    if (model->tracks.empty())
    {
        status = "Free";
        stage = model->states_labels[cell];
        return;
    }

    const size_t nk = model->nh_nstatus();
    stage = model->tracks[static_cast< size_t >(d)].nh->get_stages()[cell / nk];
    status = model->st ?
        model->st->get_stages()[cell % nk] :
        std::string("Free");

}

inline void DataBase::split_states(
    int d,
    const std::vector< epiworld_fast_uint > & cells,
    std::vector< std::string > * status,
    std::vector< std::string > * stage
) const
{

    std::string k, s;
    for (auto c : cells)
    {
        split_state(d, c, k, s);
        if (status != nullptr)
            status->push_back(k);
        if (stage != nullptr)
            stage->push_back(s);
    }

}

inline int DataBase::today_total_of(
    const std::string & name,
    bool by_stage,
    int disease
) const
{

    // A status is counted in the first disease's cells (it is the same in
    // every disease); a stage, in each disease asked for that has it
    const auto ds = by_stage ? diseases_of(disease) : std::vector< int >{0};

    int n_found = 0;
    int total = 0;
    std::string k, s;
    for (int d : ds)
    {
        const auto & r = recs[static_cast< size_t >(d)];
        bool found = false;
        for (epiworld_fast_uint i = 0u; i < r.today_total.size(); ++i)
        {
            split_state(d, i, k, s);
            if ((by_stage ? s : k) == name)
            {
                found = true;
                total += r.today_total[i];
            }
        }
        n_found += found ? 1 : 0;
    }

    if (n_found == 0)
        throw std::range_error(
            "'" + name + "' is not a " + (by_stage ? "stage" : "status") +
            " of the model."
        );

    if (n_found > 1)
        throw std::logic_error(
            "Several diseases have a stage named '" + name + "': pass the "
            "disease (disease_id())."
        );

    return total;

}

inline int DataBase::get_today_total_stage(
    const std::string & stage,
    int disease
) const
{
    return today_total_of(stage, true, disease);
}

inline int DataBase::get_today_total_status(
    const std::string & status
) const
{
    return today_total_of(status, false, all_diseases);
}

inline void DataBase::get_today_total(
    std::vector< std::string > * status,
    std::vector< std::string > * stage,
    std::vector< int > * counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    if (status != nullptr)
        status->clear();
    if (stage != nullptr)
        stage->clear();
    if (counts != nullptr)
        counts->clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        std::vector< epiworld_fast_uint > cells(r.today_total.size());
        std::iota(cells.begin(), cells.end(), 0u);
        split_states(d, cells, status, stage);
        if (counts != nullptr)
            epi_append(*counts, r.today_total);
        disease_rows(disease_names, cells.size(), d);
    }

}

inline void DataBase::get_today_virus(
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & id,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    status.clear();
    stage.clear();
    id.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    // One row per virus (of the disease) and cell
    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        std::vector< epiworld_fast_uint > cells;

        for (size_t v = 0u; v < r.today_virus.size(); ++v)
        {
            if ((recs.size() > 1u) && (disease_of_virus(static_cast< int >(v)) != d))
                continue;
            for (epiworld_fast_uint s = 0u; s < r.ncells; ++s)
            {
                cells.push_back(s);
                id.push_back(static_cast< int >(v));
                counts.push_back(r.today_virus[v][s]);
            }
        }

        split_states(d, cells, &status, &stage);
        disease_rows(disease_names, cells.size(), d);

    }

}

inline void DataBase::get_hist_total(
    std::vector< int > * date,
    std::vector< std::string > * status,
    std::vector< std::string > * stage,
    std::vector< int > * counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    if (date != nullptr)
        date->clear();
    if (status != nullptr)
        status->clear();
    if (stage != nullptr)
        stage->clear();
    if (counts != nullptr)
        counts->clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        if (date != nullptr)
            epi_append(*date, r.hist_total_date);
        split_states(d, r.hist_total_state, status, stage);
        if (counts != nullptr)
            epi_append(*counts, r.hist_total_counts);
        disease_rows(disease_names, r.hist_total_counts.size(), d);
    }

}

inline void DataBase::get_hist_virus(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    id.clear();
    status.clear();
    stage.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        epi_append(date, r.hist_virus_date);
        epi_append(id, r.hist_virus_id);
        split_states(d, r.hist_virus_state, &status, &stage);
        epi_append(counts, r.hist_virus_counts);
        disease_rows(disease_names, r.hist_virus_counts.size(), d);
    }

}

inline void DataBase::get_hist_tool(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & status,
    std::vector< std::string > & stage,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    id.clear();
    status.clear();
    stage.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        epi_append(date, r.hist_tool_date);
        epi_append(id, r.hist_tool_id);
        split_states(d, r.hist_tool_state, &status, &stage);
        epi_append(counts, r.hist_tool_counts);
        disease_rows(disease_names, r.hist_tool_counts.size(), d);
    }

}

inline void DataBase::get_hist_transition_matrix(
    std::vector< std::string > & status_from,
    std::vector< std::string > & stage_from,
    std::vector< std::string > & status_to,
    std::vector< std::string > & stage_to,
    std::vector< int > & date,
    std::vector< int > & counts,
    bool skip_zeros,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    status_from.clear();
    stage_from.clear();
    status_to.clear();
    stage_to.clear();
    date.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        const size_t n = r.hist_transition_matrix.size();
        if (n == 0u)
            continue;

        std::vector< epiworld_fast_uint > from, to;
        from.reserve(n);
        to.reserve(n);

        const size_t n_states = r.ncells;
        const size_t n_steps  = model->get_ndays();
        size_t rows = 0u;

        for (size_t step = 0u; step <= n_steps; ++step) // The final step counts
            for (size_t j = 0u; j < n_states; ++j) // Column major storage
                for (size_t i = 0u; i < n_states; ++i)
                {
                    // Retrieving the value of the day
                    int v = r.hist_transition_matrix[
                        step * n_states * n_states + // Day of the data
                        j * n_states +               // Column (to)
                        i                            // Row (from)
                        ];

                    // If we are skipping the zeros and it is zero, then don't save
                    if (skip_zeros && v == 0)
                        continue;

                    from.push_back(static_cast< epiworld_fast_uint >(i));
                    to.push_back(static_cast< epiworld_fast_uint >(j));
                    date.push_back(r.hist_total_date[step * n_states]);
                    counts.push_back(v);
                    ++rows;

                }

        split_states(d, from, &status_from, &stage_from);
        split_states(d, to, &status_to, &stage_to);
        disease_rows(disease_names, rows, d);

    }

}

inline void DataBase::get_active_cases(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & count,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    virus_id.clear();
    count.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    // Extracting useful sizes
    const size_t n_days    = model->get_ndays() + 1u; // Including day 0
    const size_t n_viruses = model->get_n_viruses();

    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];

        // Every day of each of the disease's viruses
        std::vector< int > date_d(n_days * n_viruses, 0);
        std::vector< int > id_d(n_days * n_viruses, 0);
        std::vector< int > count_d(n_days * n_viruses, 0);

        for (size_t i = 0u; i < r.hist_virus_id.size(); ++i)
        {

            // With more viruses, there are more days
            auto location = r.hist_virus_date[i] + r.hist_virus_id[i] * n_days;

            date_d[location] = r.hist_virus_date[i];
            id_d[location] = r.hist_virus_id[i];
            count_d[location] += r.hist_virus_counts[i];

        }

        size_t rows = 0u;
        for (size_t v = 0u; v < n_viruses; ++v)
        {
            if ((recs.size() > 1u) && (disease_of_virus(static_cast< int >(v)) != d))
                continue;
            for (size_t k = 0u; k < n_days; ++k)
            {
                date.push_back(date_d[k + v * n_days]);
                virus_id.push_back(id_d[k + v * n_days]);
                count.push_back(count_d[k + v * n_days]);
                ++rows;
            }
        }

        disease_rows(disease_names, rows, d);

    }

}

inline void DataBase::get_outbreak_size(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & outbreak_size,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    virus_id.clear();
    outbreak_size.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    // Extracting useful sizes
    const size_t n_days    = model->get_ndays() + 1u; // Including day 0
    const size_t n_viruses = model->get_n_viruses();

    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        std::vector< int > size_d(n_days * n_viruses, 0);

        // With more viruses, there are more days;
        for (size_t i = 0u; i < r.transmission_date.size(); ++i)
            size_d[r.transmission_date[i] + r.transmission_virus[i] * n_days] += 1;

        size_t rows = 0u;
        for (size_t v = 0u; v < n_viruses; ++v)
        {
            if ((recs.size() > 1u) && (disease_of_virus(static_cast< int >(v)) != d))
                continue;

            // The cumulative sum
            for (size_t k = 0u; k < n_days; ++k)
            {
                if (k > 0u)
                    size_d[k + v * n_days] += size_d[k - 1u + v * n_days];
                date.push_back(static_cast< int >(k + v * n_days));
                virus_id.push_back(static_cast< int >(v));
                outbreak_size.push_back(size_d[k + v * n_days]);
                ++rows;
            }
        }

        disease_rows(disease_names, rows, d);

    }

}

inline void DataBase::get_transmissions(
    std::vector<int> & date,
    std::vector<int> & source,
    std::vector<int> & target,
    std::vector<int> & virus,
    std::vector<int> & source_exposure_date,
    std::vector<int> & layer,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    source.clear();
    target.clear();
    virus.clear();
    source_exposure_date.clear();
    layer.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        epi_append(date, r.transmission_date);
        epi_append(source, r.transmission_source);
        epi_append(target, r.transmission_target);
        epi_append(virus, r.transmission_virus);
        epi_append(source_exposure_date, r.transmission_source_exposure_date);
        epi_append(layer, r.transmission_layer);
        disease_rows(disease_names, r.transmission_date.size(), d);
    }

}

inline MapVec_type<int,int> DataBase::get_reproductive_number(
    int disease
) const {

    // The keys hold the virus, so the diseases' keys never meet
    MapVec_type<int,int> map;

    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        for (size_t i = 0u; i < r.transmission_date.size(); ++i)
        {
            // Fabricating id
            std::vector< int > h = {
                r.transmission_virus[i],
                r.transmission_source[i],
                r.transmission_source_exposure_date[i]
            };

            // Adding to counter
            if (map.find(h) == map.end())
                map[h] = 1;
            else
                map[h]++;

            // The target is added
            std::vector< int > h_target = {
                r.transmission_virus[i],
                r.transmission_target[i],
                r.transmission_date[i]
            };

            map[h_target] = 0;

        }

    }

    return map;

}

inline TransitionMatrix DataBase::get_transition_probability(
    Margin margin,
    int disease,
    bool normalize
) const
{

    check_disease(disease);

    // A status matrix is the same in every disease (the first one's);
    // a pair or a stage matrix is one disease's
    int d = 0;
    if ((margin != Margin::status) && !model->tracks.empty())
        d = model->nh_disease(disease, "get_transition_probability()").id;
    else if (disease != all_diseases)
        d = disease;
    const auto & r = recs[static_cast< size_t >(d)];

    TransitionMatrix tm;
    tm.margin = margin;
    tm.normalized = normalize;
    if ((margin != Margin::status) && !model->tracks.empty())
        tm.disease = model->tracks[static_cast< size_t >(d)].nh->get_name();

    // The labels of each cell, and the margin's labels
    const size_t ns = r.ncells;
    std::vector< epiworld_fast_uint > cells(ns);
    std::iota(cells.begin(), cells.end(), 0u);
    std::vector< std::string > status, stage;
    split_states(d, cells, &status, &stage);

    std::vector< size_t > row(ns);
    if (margin == Margin::pair)
    {
        tm.status = status;
        tm.stage = stage;
        std::iota(row.begin(), row.end(), 0u);
    }
    else
    {
        const auto & key = (margin == Margin::stage) ? stage : status;
        auto & labels = (margin == Margin::stage) ? tm.stage : tm.status;
        for (size_t c = 0u; c < ns; ++c)
        {
            auto it = std::find(labels.begin(), labels.end(), key[c]);
            row[c] = static_cast< size_t >(it - labels.begin());
            if (it == labels.end())
                labels.push_back(key[c]);
        }
    }

    // The counts over the run's days (the record after the last day is not
    // a day of transitions), summed over the other coordinate
    const size_t n = tm.size();
    tm.values.assign(n * n, 0.0);
    const size_t recorded = (ns == 0u) ? 0u : r.hist_transition_matrix.size() / (ns * ns);
    const size_t ndays = std::min(static_cast< size_t >(model->get_ndays()), recorded);
    for (size_t t = 0u; t < ndays; ++t)
        for (size_t to = 0u; to < ns; ++to)
            for (size_t from = 0u; from < ns; ++from)
                tm.values[row[from] + row[to] * n] += static_cast< epiworld_double >(
                    r.hist_transition_matrix[t * ns * ns + to * ns + from]
                );

    if (normalize)
        for (size_t from = 0u; from < n; ++from)
        {
            epiworld_double total = 0.0;
            for (size_t to = 0u; to < n; ++to)
                total += tm.values[from + to * n];
            if (total == 0.0)
                continue;
            for (size_t to = 0u; to < n; ++to)
                tm.values[from + to * n] /= total;
        }

    return tm;

}

inline std::string TransitionMatrix::label(size_t i) const
{
    if (margin == Margin::pair)
        return status.at(i) + " / " + stage.at(i);
    return (margin == Margin::stage) ? stage.at(i) : status.at(i);
}

inline void TransitionMatrix::print() const
{

    const size_t n = size();

    size_t nchar = 0u;
    for (size_t i = 0u; i < n; ++i)
        nchar = std::max(nchar, label(i).length());

    std::string fmt = " - %-" + std::to_string(nchar) + "s";
    std::string fmt_entry = " % 4.2f";
    if (!normalized)
    {
        size_t width = 0u;
        for (auto v : values)
            width = std::max(width, std::to_string(static_cast< long long >(v)).length());
        fmt_entry = " % " + std::to_string(width) + ".0f";
    }

    printf_epiworld("\nTransition Probabilities:\n");
    for (size_t from = 0u; from < n; ++from)
    {
        printf_epiworld(fmt.c_str(), label(from).c_str());
        for (size_t to = 0u; to < n; ++to)
        {
            const epiworld_double v = values[from + to * n];
            // (printf_epiworld is two statements: braces)
            if (std::isnan(v) || (v < 1e-10)) {
                printf_epiworld("     -");
            } else {
                printf_epiworld(fmt_entry.c_str(), v);
            }
        }
        printf_epiworld("\n");
    }

    printf_epiworld("\n");

}

inline void DataBase::get_generation_time(
    std::vector< int > & agent_id,
    std::vector< int > & virus_id,
    std::vector< int > & time,
    std::vector< int > & gentime,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    agent_id.clear();
    virus_id.clear();
    time.clear();
    gentime.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        const size_t nevents = r.transmission_date.size();

        // Iterating through the individuals
        for (size_t i = 0u; i < nevents; ++i)
        {
            const int agent_id_i = r.transmission_target[i];
            agent_id.push_back(agent_id_i);
            virus_id.push_back(r.transmission_virus[i]);
            time.push_back(r.transmission_date[i]);

            bool found = false;
            for (size_t j = i; j < nevents; ++j)
            {

                if (r.transmission_source[j] == agent_id_i)
                {
                    gentime.push_back(r.transmission_date[j] - r.transmission_date[i]);
                    found = true;
                    break;
                }

            }

            // If there's no transmission, we set the generation time to
            // minus 1;
            if (!found)
                gentime.push_back(-1);

        }

        disease_rows(disease_names, nevents, d);

    }

}

inline void DataBase::get_hospitalizations(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & tool_id,
    std::vector<int> & count,
    std::vector<double> & weight,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    date.clear();
    virus_id.clear();
    tool_id.clear();
    count.clear();
    weight.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    // Get the number of days from the model and add 1 since days are 0-indexed
    // today() returns the last day, so we need today() + 1 for the full count
    const int ndays = model->today() + 1;
    for (int d : diseases_of(disease))
    {
        std::vector< int > date_d, virus_d, tool_d, count_d;
        std::vector< double > weight_d;
        recs[static_cast< size_t >(d)].m_hospitalizations.get(
            ndays, date_d, virus_d, tool_d, count_d, weight_d
        );
        epi_append(date, date_d);
        epi_append(virus_id, virus_d);
        epi_append(tool_id, tool_d);
        epi_append(count, count_d);
        epi_append(weight, weight_d);
        disease_rows(disease_names, date_d.size(), d);
    }

}

inline void DataBase::get_hist_virus_transitions(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & transition,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    if (model->tracks.empty())
        throw std::logic_error(
            "get_hist_virus_transitions() needs a model with a natural history."
        );

    date.clear();
    id.clear();
    transition.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    for (int d : diseases_of(disease))
    {
        const auto & r = recs[static_cast< size_t >(d)];
        const auto & nh = *model->tracks[static_cast< size_t >(d)].nh;
        epi_append(date, r.hist_transition_virus_date);
        epi_append(id, r.hist_transition_virus_id);
        epi_append(counts, r.hist_transition_virus_counts);
        for (int t : r.hist_transition_virus_transition)
            transition.push_back(nh.get_transition(t).name);
        disease_rows(disease_names, r.hist_transition_virus_id.size(), d);
    }

}

inline void DataBase::get_hist_transitions(
    std::vector< int > & date,
    std::vector< std::string > & transition,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    if (model->tracks.empty())
        throw std::logic_error(
            "get_hist_transitions() needs a model with a natural history."
        );

    date.clear();
    transition.clear();
    counts.clear();
    if (disease_names != nullptr)
        disease_names->clear();

    // Each disease's transitions; the statuses', recorded once after the
    // first disease's, belong to every disease: with several diseases they
    // come with the disease "" (and only when every disease is asked for)
    const bool several = recs.size() > 1u;
    for (int d : diseases_of(disease))
    {

        const auto & r = recs[static_cast< size_t >(d)];
        const auto & nh = *model->tracks[static_cast< size_t >(d)].nh;
        const int n_stage = static_cast< int >(nh.get_n_transitions());

        size_t rows = 0u;
        size_t status_rows = 0u;
        for (size_t k = 0u; k < r.hist_transition_id.size(); ++k)
        {
            const int t = r.hist_transition_id[k];
            if ((t >= n_stage) && several)
                continue;
            date.push_back(r.hist_transition_date[k]);
            counts.push_back(r.hist_transition_counts[k]);
            transition.push_back(
                (t < n_stage) ? nh.get_transition(t).name :
                    model->st->get_transition(t - n_stage).name
            );
            ++rows;
        }
        disease_rows(disease_names, rows, d);

        // Several diseases: the statuses' rows, once
        if (several && (d == 0) && (disease == all_diseases))
        {
            for (size_t k = 0u; k < r.hist_transition_id.size(); ++k)
            {
                const int t = r.hist_transition_id[k];
                if (t < n_stage)
                    continue;
                date.push_back(r.hist_transition_date[k]);
                counts.push_back(r.hist_transition_counts[k]);
                transition.push_back(model->st->get_transition(t - n_stage).name);
                ++status_rows;
            }
            disease_rows(disease_names, status_rows, -1);
        }

    }

}

inline void DataBase::check_disease(int disease) const
{

    const int n = static_cast< int >(model->tracks.size());
    if ((disease == all_diseases) || ((disease >= 0) && (disease < n)))
        return;

    std::string names;
    for (int d = 0; d < n; ++d)
        names += (d ? ", " : "") + std::to_string(d) + " (\"" +
            model->tracks[static_cast< size_t >(d)].nh->get_name() + "\")";

    throw std::range_error(
        "There is no disease " + std::to_string(disease) + " in the model. "
        "Its diseases are: " + (n ? names : std::string("none")) +
        "; all_diseases (-1) asks for every disease."
    );

}

inline void DataBase::disease_rows(
    std::vector< std::string > * names,
    size_t n,
    int d
) const
{

    if (names == nullptr)
        return;

    names->insert(
        names->end(), n,
        ((d < 0) || model->tracks.empty()) ? std::string() :
            model->tracks[static_cast< size_t >(d)].nh->get_name()
    );

}

#endif
