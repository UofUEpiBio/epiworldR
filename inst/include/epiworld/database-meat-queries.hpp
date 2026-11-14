#ifndef EPIWORLD_DATABASE_MEAT_QUERIES_HPP
#define EPIWORLD_DATABASE_MEAT_QUERIES_HPP

/**
 * @file database-meat-queries.hpp
 * @brief Reading what the database recorded: counts, histories,
 * transitions, transmissions, reproductive numbers and hospitalizations.
 */

inline void DataBase::split_state(
    epiworld_fast_uint state,
    std::string & status,
    std::string & stage
) const
{

    if (model->tracks.empty())
    {
        status = "Free";
        stage = model->states_labels[state];
        return;
    }

    const unsigned int c = static_cast< unsigned int >(state);
    stage = model->tracks[0u].nh->get_stages()[static_cast< size_t >(model->nh_stage_of(c))];
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
    for (epiworld_fast_uint i = 0u; i < rec.today_total.size(); ++i)
    {
        split_state(i, k, s);
        if ((by_stage ? s : k) == name)
        {
            found = true;
            total += rec.today_total[i];
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
    const std::string & stage,
    int disease
) const
{
    check_disease(disease);
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
    std::vector< int > * counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    check_disease(disease);


    std::vector< epiworld_fast_uint > states(rec.today_total.size());
    for (size_t i = 0u; i < states.size(); ++i)
        states[i] = static_cast< epiworld_fast_uint >(i);

    split_states(states, status, stage);

    if (counts != nullptr)
        *counts = rec.today_total;

    disease_rows(disease_names, rec.today_total.size());

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

    check_disease(disease);


    // One row per virus and state
    const size_t n_viruses = rec.today_virus.size();
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
            counts.push_back(rec.today_virus[v][s]);
        }

    split_states(states, &status, &stage);

    disease_rows(disease_names, id.size());

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

    check_disease(disease);


    if (date != nullptr)
        *date = rec.hist_total_date;

    split_states(rec.hist_total_state, status, stage);

    if (counts != nullptr)
        *counts = rec.hist_total_counts;

    disease_rows(disease_names, rec.hist_total_counts.size());

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

    check_disease(disease);


    date = rec.hist_virus_date;
    id = rec.hist_virus_id;
    split_states(rec.hist_virus_state, &status, &stage);
    counts = rec.hist_virus_counts;

    disease_rows(disease_names, counts.size());

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

    check_disease(disease);


    date = rec.hist_tool_date;
    id = rec.hist_tool_id;
    split_states(rec.hist_tool_state, &status, &stage);
    counts = rec.hist_tool_counts;

    disease_rows(disease_names, counts.size());

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

    check_disease(disease);


    size_t n = rec.hist_transition_matrix.size();

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
                int v = rec.hist_transition_matrix[
                    step * n_states * n_states + // Day of the data
                    j * n_states +               // Column (to)
                    i                            // Row (from)
                    ];

                // If we are skipping the zeros and it is zero, then don't save
                if (skip_zeros && v == 0)
                    continue;
                                
                from.push_back(static_cast< epiworld_fast_uint >(i));
                to.push_back(static_cast< epiworld_fast_uint >(j));
                date.push_back(rec.hist_total_date[step * n_states]);
                counts.push_back(v);

            }

        }
    }

    split_states(from, &status_from, &stage_from);
    split_states(to, &status_to, &stage_to);

    disease_rows(disease_names, counts.size());

}

inline void DataBase::get_active_cases(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & count,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    check_disease(disease);


    // Extracting useful sizes
    size_t n_days    = model->get_ndays() + 1u; // Including day 0
    size_t n_viruses = model->get_n_viruses();
    
    // Making room
    date.assign(n_days * n_viruses, 0);
    virus_id.assign(n_days * n_viruses, 0);
    count.assign(n_days * n_viruses, 0);

    for (size_t i = 0u; i < rec.hist_virus_id.size(); ++i)
    {

        // With more viruses, there are more days
        auto location = rec.hist_virus_date[i] + rec.hist_virus_id[i] * n_days;

        date[location] = rec.hist_virus_date[i];
        virus_id[location] = rec.hist_virus_id[i];
        count[location] += rec.hist_virus_counts[i];

    }

    disease_rows(disease_names, count.size());

}

inline void DataBase::get_outbreak_size(
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & outbreak_size,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    check_disease(disease);


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
    for (size_t i = 0u; i < rec.transmission_date.size(); ++i)
    {
        outbreak_size[
            rec.transmission_date[i] + rec.transmission_virus[i] * n_days
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

    disease_rows(disease_names, outbreak_size.size());

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

    check_disease(disease);


    size_t nevents = rec.transmission_date.size();

    date.resize(nevents);
    source.resize(nevents);
    target.resize(nevents);
    virus.resize(nevents);
    source_exposure_date.resize(nevents);
    layer.resize(nevents);

    get_transmissions_into(
        date.data(),
        source.data(),
        target.data(),
        virus.data(),
        source_exposure_date.data(),
        layer.data()
    );

    disease_rows(disease_names, date.size());

}

inline void DataBase::get_transmissions_into(
    int * date,
    int * source,
    int * target,
    int * virus,
    int * source_exposure_date,
    int * layer
) const 
{

    size_t nevents = rec.transmission_date.size();

    for (size_t i = 0u; i < nevents; ++i)
    {

        *(date + i) = rec.transmission_date.at(i);
        *(source + i) = rec.transmission_source.at(i);
        *(target + i) = rec.transmission_target.at(i);
        *(virus + i) = rec.transmission_virus.at(i);
        *(source_exposure_date + i) = rec.transmission_source_exposure_date.at(i);
        *(layer + i) = rec.transmission_layer.at(i);

    }

}

inline MapVec_type<int,int> DataBase::get_reproductive_number(
    int disease
) const {

    check_disease(disease);

    // Checking size
    MapVec_type<int,int> map;

    // Number of digits of maxid
    for (size_t i = 0u; i < rec.transmission_date.size(); ++i)
    {
        // Fabricating id
        std::vector< int > h = {
            rec.transmission_virus[i],
            rec.transmission_source[i],
            rec.transmission_source_exposure_date[i]
        };

        // Adding to counter
        if (map.find(h) == map.end())
            map[h] = 1;
        else
            map[h]++;

        // The target is added
        std::vector< int > h_target = {
            rec.transmission_virus[i],
            rec.transmission_target[i],
            rec.transmission_date[i]
        };

        map[h_target] = 0;
        
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

    TransitionMatrix tm;
    tm.margin = margin;
    tm.normalized = normalize;
    if ((margin != Margin::status) && !model->tracks.empty())
        tm.disease = model->tracks[0u].nh->get_name();

    // The labels of each state (cell), and the margin's labels
    const size_t ns = model->nstates;
    std::vector< epiworld_fast_uint > cells(ns);
    std::iota(cells.begin(), cells.end(), 0u);
    std::vector< std::string > status, stage;
    split_states(cells, &status, &stage);

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
    const size_t recorded = (ns == 0u) ? 0u : rec.hist_transition_matrix.size() / (ns * ns);
    const size_t ndays = std::min(static_cast< size_t >(model->get_ndays()), recorded);
    for (size_t t = 0u; t < ndays; ++t)
        for (size_t to = 0u; to < ns; ++to)
            for (size_t from = 0u; from < ns; ++from)
                tm.values[row[from] + row[to] * n] += static_cast< epiworld_double >(
                    rec.hist_transition_matrix[t * ns * ns + to * ns + from]
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

    check_disease(disease);

    
    size_t nevents = rec.transmission_date.size();

    agent_id.reserve(nevents);
    virus_id.reserve(nevents);
    time.reserve(nevents);
    gentime.reserve(nevents);

    // Iterating through the individuals
    for (size_t i = 0u; i < nevents; ++i)
    {
        int agent_id_i = rec.transmission_target[i];
        agent_id.push_back(agent_id_i);
        virus_id.push_back(rec.transmission_virus[i]);
        time.push_back(rec.transmission_date[i]);

        bool found = false;
        for (size_t j = i; j < nevents; ++j)
        {

            if (rec.transmission_source[j] == agent_id_i)
            {
                gentime.push_back(rec.transmission_date[j] - time[i]);
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

    disease_rows(disease_names, agent_id.size());

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

    check_disease(disease);

    // Get the number of days from the model and add 1 since days are 0-indexed
    // today() returns the last day, so we need today() + 1 for the full count
    int ndays = model->today() + 1;
    rec.m_hospitalizations.get(ndays, date, virus_id, tool_id, count, weight);

    disease_rows(disease_names, date.size());

}

inline void DataBase::get_hist_transitions(
    std::vector< int > & date,
    std::vector< std::string > & transition,
    std::vector< int > & counts,
    int disease,
    std::vector< std::string > * disease_names
) const
{

    check_disease(disease);


    if (model->tracks.empty())
        throw std::logic_error(
            "get_hist_transitions() needs a model with a natural history."
        );

    date = rec.hist_transition_date;
    counts = rec.hist_transition_counts;
    // Stage transitions first, then the status ones
    const int n_stage = static_cast< int >(model->tracks[0u].nh->get_n_transitions());
    transition.resize(rec.hist_transition_id.size());
    for (size_t k = 0u; k < rec.hist_transition_id.size(); ++k)
    {
        const int t = rec.hist_transition_id[k];
        transition[k] = (t < n_stage) ?
            model->tracks[0u].nh->get_transition(t).name :
            model->st->get_transition(t - n_stage).name;
    }

    disease_rows(disease_names, counts.size());

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
    size_t n
) const
{

    if (names == nullptr)
        return;

    // One disease for now (several: NH-5b-2, several-diseases.md)
    names->assign(
        n,
        model->tracks.empty() ? std::string() : model->tracks[0u].nh->get_name()
    );

}

#endif
