#ifndef EPIWORLD_DATABASE_MEAT_HPP
#define EPIWORLD_DATABASE_MEAT_HPP

template<typename TSeq>
inline void DataBase<TSeq>::reset()
{

    // Initializing the counts
    today_total.resize(model->nstates);
    std::fill(today_total.begin(), today_total.end(), 0);
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

    
    transition_matrix.resize(model->nstates * model->nstates);
    std::fill(transition_matrix.begin(), transition_matrix.end(), 0);
    for (size_t s = 0u; s < model->nstates; ++s)
        transition_matrix[s + s * model->nstates] = today_total[s];

    hist_virus_date.clear();
    hist_virus_id.clear();
    hist_virus_state.clear();
    hist_virus_counts.clear();

    hist_tool_date.clear();
    hist_tool_id.clear();
    hist_tool_state.clear();
    hist_tool_counts.clear();    

    today_virus.resize(get_n_viruses());
    std::fill(today_virus.begin(), today_virus.begin(), std::vector<int>(model->nstates, 0));

    today_tool.resize(get_n_tools());
    std::fill(today_tool.begin(), today_tool.begin(), std::vector<int>(model->nstates, 0));

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

    return;

}

template<typename TSeq>
inline DataBase<TSeq>::DataBase(const DataBase<TSeq> & db) :
    virus_id(db.virus_id),
    virus_name(db.virus_name),
    virus_sequence(db.virus_sequence),
    virus_origin_date(db.virus_origin_date),
    virus_parent_id(db.virus_parent_id),
    tool_id(db.tool_id),
    tool_name(db.tool_name),
    tool_sequence(db.tool_sequence),
    tool_origin_date(db.tool_origin_date),
    seq_hasher(db.seq_hasher),
    seq_writer(db.seq_writer),
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
    transition_matrix(db.transition_matrix),
    user_data(nullptr)
{}

// DataBase<TSeq> & DataBase<TSeq>::operator=(const DataBase<TSeq> & m)
// {

// }

template<typename TSeq>
inline Model<TSeq> * DataBase<TSeq>::get_model() {
    return model;
}

template<typename TSeq>
inline const std::vector< TSeq > & DataBase<TSeq>::get_sequence() const {
    return virus_sequence;
}

template<typename TSeq>
inline void DataBase<TSeq>::record() 
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
        for (auto & p : virus_id)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                hist_virus_date.push_back(model->today());
                hist_virus_id.push_back(p.second);
                hist_virus_state.push_back(s);
                hist_virus_counts.push_back(today_virus[p.second][s]);

            }

        }

        // Recording tool's history
        for (auto & p : tool_id)
        {

            for (epiworld_fast_uint s = 0u; s < model->nstates; ++s)
            {

                hist_tool_date.push_back(model->today());
                hist_tool_id.push_back(p.second);
                hist_tool_state.push_back(s);
                hist_tool_counts.push_back(today_tool[p.second][s]);

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

        for (auto cell : transition_matrix)
            hist_transition_matrix.push_back(cell);

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

template<typename TSeq>
inline void DataBase<TSeq>::record_virus(Virus<TSeq> & v)
{

    // If no sequence, then need to add one. This is regardless of the case
    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        if (v.get_sequence() == -1)
            v.set_sequence(default_sequence<TSeq>(
                static_cast<int>(virus_name.size())
                ));
    }
    else
    {
        if (v.get_sequence() == nullptr)
            v.set_sequence(default_sequence<TSeq>(
                static_cast<int>(virus_name.size())
                ));        
    }

    // Negative id -> virus hasn't been recorded
    if (v.get_id() < 0)
    {

        epiworld_fast_uint new_id = virus_id.size();
        virus_name.push_back(v.get_name());

        // Generating the hash
        std::vector< int > hash;
        EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
        {
            hash = seq_hasher(v.get_sequence());
            virus_id[hash] = new_id;
            virus_sequence.push_back(v.get_sequence());
        }
        else
        {
            hash = seq_hasher(*v.get_sequence());
            virus_id[hash] = new_id;
            virus_sequence.push_back(*v.get_sequence());
        }


        virus_origin_date.push_back(model->today());
        
        virus_parent_id.push_back(v.get_id()); // Must be -99
        
        today_virus.push_back({});
        today_virus[new_id].resize(model->nstates, 0);
       
        // Updating the variant
        v.set_id(new_id);
        v.set_date(model->today());

        today_total_nviruses_active++;

    }
    else
    { // In this case, the virus is already on record, need to make sure
             // The new sequence is new.

        // Updating registry
        std::vector< int > hash;
        EPI_IF_TSEQ_LESS_EQ_INT(TSeq)
        {
            hash = seq_hasher(v.get_sequence());
        }
        else
        {
            hash = seq_hasher(*v.get_sequence());
        }
        epiworld_fast_uint old_id = v.get_id();
        epiworld_fast_uint new_id;

        // If the sequence is new, then it means that the
        if (virus_id.find(hash) == virus_id.end())
        {

            new_id = virus_id.size();
            virus_id[hash] = new_id;
            virus_name.push_back(v.get_name());

            EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
            {
                virus_sequence.push_back(v.get_sequence());
            }
            else
            {
                virus_sequence.push_back(*v.get_sequence());
            }

            virus_origin_date.push_back(model->today());
            
            virus_parent_id.push_back(old_id);
            
            today_virus.push_back({});
            today_virus[new_id].resize(model->nstates, 0);
        
            // Updating the variant
            v.set_id(new_id);
            v.set_date(model->today());

            today_total_nviruses_active++;

        } else {

            // Finding the id
            new_id = virus_id[hash];

            // Reflecting the change
            v.set_id(new_id);
            v.set_date(virus_origin_date[new_id]);

        }

        // Moving statistics (only if we are affecting an individual)
        if (v.get_agent() != nullptr)
        {
            // Correcting math
            epiworld_fast_uint tmp_state = v.get_agent()->get_state();
            today_virus[old_id][tmp_state]--;
            today_virus[new_id][tmp_state]++;

        }

    }
    
    return;

} 

template<typename TSeq>
inline void DataBase<TSeq>::record_tool(Tool<TSeq> & t)
{

    EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
    {
        if (t.get_sequence() == -1)
            t.set_sequence(default_sequence<TSeq>(
                static_cast<int>(tool_name.size())
            ));
    }
    else
    {
        if (t.get_sequence() == nullptr)
            t.set_sequence(default_sequence<TSeq>(
                static_cast<int>(tool_name.size())
            ));
    }

    if (t.get_id() < 0) 
    {

        epiworld_fast_uint new_id = tool_id.size();
        tool_name.push_back(t.get_name());

        std::vector< int > hash;
        EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
        {
            hash = seq_hasher(t.get_sequence());
            tool_id[hash] = new_id;
            tool_sequence.push_back(t.get_sequence());
        }
        else
        {
            hash = seq_hasher(*t.get_sequence());
            tool_id[hash] = new_id;
            tool_sequence.push_back(*t.get_sequence());

        }
        tool_origin_date.push_back(model->today());
                
        today_tool.push_back({});
        today_tool[new_id].resize(model->nstates, 0);

        // Updating the tool
        t.set_id(new_id);
        t.set_date(model->today());

    } else {

        // Updating registry
        std::vector< int > hash;
        EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
        {
            hash = seq_hasher(t.get_sequence());
        }
        else
        {
            hash = seq_hasher(*t.get_sequence());
        }
        epiworld_fast_uint old_id = t.get_id();
        epiworld_fast_uint new_id;
        
        if (tool_id.find(hash) == tool_id.end())
        {

            new_id = tool_id.size();
            tool_id[hash] = new_id;
            tool_name.push_back(t.get_name());

            EPI_IF_TSEQ_LESS_EQ_INT( TSeq )
            {
                tool_sequence.push_back(t.get_sequence());
            }
            else
            {
                tool_sequence.push_back(*t.get_sequence());
            }
            
            tool_origin_date.push_back(model->today());
                    
            today_tool.push_back({});
            today_tool[new_id].resize(model->nstates, 0);

            // Updating the tool
            t.set_id(new_id);
            t.set_date(model->today());

        } else {

            // Finding the id
            new_id = tool_id[hash];

            // Reflecting the change
            t.set_id(new_id);
            t.set_date(tool_origin_date[new_id]);

        }

        // Moving statistics (only if we are affecting an individual)
        if (t.get_agent() != nullptr)
        {
            // Correcting math
            epiworld_fast_uint tmp_state = t.get_agent()->get_state();
            today_tool[old_id][tmp_state]--;
            today_tool[new_id][tmp_state]++;

        }

    }

    
    
    return;
} 

template<typename TSeq>
inline size_t DataBase<TSeq>::size() const
{
    return virus_id.size();
}

template<typename TSeq>
inline void DataBase<TSeq>::update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo
) {

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

template<typename TSeq>
inline void DataBase<TSeq>::update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
) {

    today_virus[virus_id][prev_state]--;
    today_virus[virus_id][new_state]++;

    return;
    
}

template<typename TSeq>
inline void DataBase<TSeq>::update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
) {


    today_tool[tool_id][prev_state]--;    
    today_tool[tool_id][new_state]++;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::record_transition(
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

template<typename TSeq>
inline int DataBase<TSeq>::get_today_total(
    std::string what
) const
{

    for (auto i = 0u; i < model->states_labels.size(); ++i)
    {
        if (model->states_labels[i] == what)
            return today_total[i];
    }

    throw std::range_error("The value '" + what + "' is not in the model.");

}

template<typename TSeq>
inline void DataBase<TSeq>::get_today_total(
    std::vector< std::string > * state,
    std::vector< int > * counts
) const
{
    if (state != nullptr)
        (*state) = model->states_labels;

    if (counts != nullptr)
        *counts = today_total;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_today_virus(
    std::vector< std::string > & state,
    std::vector< int > & id,
    std::vector< int > & counts
    ) const
{
      
    state.resize(today_virus.size(), "");
    id.resize(today_virus.size(), 0);
    counts.resize(today_virus.size(),0);

    int n = 0u;
    for (epiworld_fast_uint v = 0u; v < today_virus.size(); ++v)
        for (epiworld_fast_uint s = 0u; s < model->states_labels.size(); ++s)
        {
            state[n]   = model->states_labels[s];
            id[n]       = static_cast<int>(v);
            counts[n++] = today_virus[v][s];

        }

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_total(
    std::vector< int > * date,
    std::vector< std::string > * state,
    std::vector< int > * counts
) const
{

    if (date != nullptr)
        *date = hist_total_date;

    if (state != nullptr)
    {
        state->resize(hist_total_state.size(), "");
        for (epiworld_fast_uint i = 0u; i < hist_total_state.size(); ++i)
            state->operator[](i) = model->states_labels[hist_total_state[i]];
    }

    if (counts != nullptr)
        *counts = hist_total_counts;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_virus(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & state,
    std::vector< int > & counts
) const {

    date = hist_virus_date;
    std::vector< std::string > labels;
    labels = model->states_labels;
    
    id = hist_virus_id;
    state.resize(hist_virus_state.size(), "");
    for (epiworld_fast_uint i = 0u; i < hist_virus_state.size(); ++i)
        state[i] = labels[hist_virus_state[i]];

    counts = hist_virus_counts;

    return;

}


template<typename TSeq>
inline void DataBase<TSeq>::get_hist_tool(
    std::vector< int > & date,
    std::vector< int > & id,
    std::vector< std::string > & state,
    std::vector< int > & counts
) const {

    date = hist_tool_date;
    std::vector< std::string > labels;
    labels = model->states_labels;
    
    id = hist_tool_id;
    state.resize(hist_tool_state.size(), "");
    for (size_t i = 0u; i < hist_tool_state.size(); ++i)
        state[i] = labels[hist_tool_state[i]];

    counts = hist_tool_counts;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_today_transition_matrix(
    std::vector< int > & counts
) const
{

    counts = transition_matrix;

    return;

}

template<typename TSeq>
inline void DataBase<TSeq>::get_hist_transition_matrix(
    std::vector< std::string > & state_from,
    std::vector< std::string > & state_to,
    std::vector< int > & date,
    std::vector< int > & counts,
    bool skip_zeros
) const
{

    size_t n = this->hist_transition_matrix.size();
    
    // Clearing the previous vectors
    state_from.clear();
    state_to.clear();
    date.clear();
    counts.clear();

    // Reserving space
    state_from.reserve(n);
    state_to.reserve(n);
    date.reserve(n);
    counts.reserve(n);

    size_t n_states = model->nstates;
    size_t n_steps  = model->get_ndays();

    // If n is zero, then we are done
    if (n == 0u)
        return;

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
                                
                state_from.push_back(model->states_labels[i]);
                state_to.push_back(model->states_labels[j]);
                date.push_back(hist_total_date[step * n_states]);
                counts.push_back(v);

            }

        }
    }

    return;


}

template<typename TSeq>
inline void DataBase<TSeq>::get_transmissions(
    std::vector<int> & date,
    std::vector<int> & source,
    std::vector<int> & target,
    std::vector<int> & virus,
    std::vector<int> & source_exposure_date
) const 
{

    size_t nevents = transmission_date.size();

    date.resize(nevents);
    source.resize(nevents);
    target.resize(nevents);
    virus.resize(nevents);
    source_exposure_date.resize(nevents);

    get_transmissions(
        &date[0u],
        &source[0u],
        &target[0u],
        &virus[0u],
        &source_exposure_date[0u]
    );

}

template<typename TSeq>
inline void DataBase<TSeq>::get_transmissions(
    int * date,
    int * source,
    int * target,
    int * virus,
    int * source_exposure_date
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

    }

}

template<typename TSeq>
inline void DataBase<TSeq>::write_data(
    std::string fn_virus_info,
    std::string fn_virus_hist,
    std::string fn_tool_info,
    std::string fn_tool_hist,
    std::string fn_total_hist,
    std::string fn_transmission,
    std::string fn_transition,
    std::string fn_reproductive_number,
    std::string fn_generation_time
) const
{

    if (fn_virus_info != "")
    {
        std::ofstream file_virus_info(fn_virus_info, std::ios_base::out);

        // Check if the file exists and throw an error if it doesn't
        if (!file_virus_info)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_virus_info +
                "\" for writing.")
                ;
        }


        file_virus_info <<
        #ifdef EPI_DEBUG
            "thread" << "virus_id " << "virus " << "virus_sequence " << "date_recorded " << "parent\n";
        #else
            "virus_id " << "virus " << "virus_sequence " << "date_recorded " << "parent\n";
        #endif

        for (const auto & v : virus_id)
        {
            int id = v.second;
            file_virus_info <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                id << " \"" <<
                virus_name[id] << "\" " <<
                seq_writer(virus_sequence[id]) << " " <<
                virus_origin_date[id] << " " <<
                virus_parent_id[id] << "\n";
        }

    }

    if (fn_virus_hist != "")
    {
        std::ofstream file_virus(fn_virus_hist, std::ios_base::out);
        
        // Repeat the same error if the file doesn't exists
        if (!file_virus)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_virus_hist +
                "\" for writing.")
                ;
        }

        file_virus <<
            #ifdef EPI_DEBUG
            "thread "<< "date " << "virus_id " << "virus " << "state " << "n\n";
            #else
            "date " << "virus_id " << "virus " << "state " << "n\n";
            #endif

        for (epiworld_fast_uint i = 0; i < hist_virus_id.size(); ++i)
            file_virus <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_virus_date[i] << " " <<
                hist_virus_id[i] << " \"" <<
                virus_name[hist_virus_id[i]] << "\" " <<
                model->states_labels[hist_virus_state[i]] << " " <<
                hist_virus_counts[i] << "\n";
    }

    if (fn_tool_info != "")
    {
        std::ofstream file_tool_info(fn_tool_info, std::ios_base::out);

        // Repeat the same error if the file doesn't exists
        if (!file_tool_info)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_tool_info +
                "\" for writing.")
                ;
        }

        file_tool_info <<
            #ifdef EPI_DEBUG
            "thread " << 
            #endif
            "id " << "tool_name " << "tool_sequence " << "date_recorded\n";

        for (const auto & t : tool_id)
        {
            int id = t.second;
            file_tool_info <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                id << " \"" <<
                tool_name[id] << "\" " <<
                seq_writer(tool_sequence[id]) << " " <<
                tool_origin_date[id] << "\n";
        }

    }

    if (fn_tool_hist != "")
    {
        std::ofstream file_tool_hist(fn_tool_hist, std::ios_base::out);

        // Repeat the same error if the file doesn't exists
        if (!file_tool_hist)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_tool_hist +
                "\" for writing.")
                ;
        }
        
        file_tool_hist <<
            #ifdef EPI_DEBUG
            "thread " << 
            #endif
            "date " << "id " << "state " << "n\n";

        for (epiworld_fast_uint i = 0; i < hist_tool_id.size(); ++i)
            file_tool_hist <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_tool_date[i] << " " <<
                hist_tool_id[i] << " " <<
                model->states_labels[hist_tool_state[i]] << " " <<
                hist_tool_counts[i] << "\n";
    }

    if (fn_total_hist != "")
    {
        std::ofstream file_total(fn_total_hist, std::ios_base::out);

        // Repeat the same error if the file doesn't exists
        if (!file_total)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_total_hist +
                "\" for writing.")
                ;
        }

        file_total <<
            #ifdef EPI_DEBUG
            "thread " << 
            #endif
            "date " << "nviruses " << "state " << "counts\n";

        for (epiworld_fast_uint i = 0; i < hist_total_date.size(); ++i)
            file_total <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                hist_total_date[i] << " " <<
                hist_total_nviruses_active[i] << " \"" <<
                model->states_labels[hist_total_state[i]] << "\" " << 
                hist_total_counts[i] << "\n";
    }

    if (fn_transmission != "")
    {
        std::ofstream file_transmission(fn_transmission, std::ios_base::out);

        // Repeat the same error if the file doesn't exists
        if (!file_transmission)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_transmission +
                "\" for writing.")
                ;
        }

        file_transmission <<
            #ifdef EPI_DEBUG
            "thread " << 
            #endif
            "date " << "virus_id virus " << "source_exposure_date " << "source " << "target\n";

        for (epiworld_fast_uint i = 0; i < transmission_target.size(); ++i)
            file_transmission <<
                #ifdef EPI_DEBUG
                EPI_GET_THREAD_ID() << " " <<
                #endif
                transmission_date[i] << " " <<
                transmission_virus[i] << " \"" <<
                virus_name[transmission_virus[i]] << "\" " <<
                transmission_source_exposure_date[i] << " " <<
                transmission_source[i] << " " <<
                transmission_target[i] << "\n";
                
    }

    if (fn_transition != "")
    {
        std::ofstream file_transition(fn_transition, std::ios_base::out);

        // Repeat the same error if the file doesn't exists
        if (!file_transition)
        {
            throw std::runtime_error(
                "Could not open file \"" + fn_transition +
                "\" for writing.")
                ;
        }

        file_transition <<
            #ifdef EPI_DEBUG
            "thread " << 
            #endif
            "date " << "from " << "to " << "counts\n";

        int ns = model->nstates;

        for (int i = 0; i <= model->today(); ++i)
        {

            // Skipping the zeros
            if (hist_transition_matrix[i * (ns * ns)] == 0)
                continue;

            for (int from = 0u; from < ns; ++from)
                for (int to = 0u; to < ns; ++to)
                    file_transition <<
                        #ifdef EPI_DEBUG
                        EPI_GET_THREAD_ID() << " " <<
                        #endif
                        i << " \"" <<
                        model->states_labels[from] << "\" \"" <<
                        model->states_labels[to] << "\" " <<
                        hist_transition_matrix[i * (ns * ns) + to * ns + from] << "\n";
                
        }
                
    }

    if (fn_reproductive_number != "")
        reproductive_number(fn_reproductive_number);

    if (fn_generation_time != "")
        generation_time(fn_generation_time);

}

template<typename TSeq>
inline void DataBase<TSeq>::record_transmission(
    int i,
    int j,
    int virus,
    int i_expo_date
) {

    transmission_date.push_back(model->today());
    transmission_source.push_back(i);
    transmission_target.push_back(j);
    transmission_virus.push_back(virus);
    transmission_source_exposure_date.push_back(i_expo_date);

}

template<typename TSeq>
inline size_t DataBase<TSeq>::get_n_viruses() const
{
    return virus_id.size();
}

template<typename TSeq>
inline size_t DataBase<TSeq>::get_n_tools() const
{
    return tool_id.size();
}


template<typename TSeq>
inline void DataBase<TSeq>::set_user_data(
    std::vector< std::string > names
)
{
    user_data = UserData<TSeq>(names);
    user_data.model = model;
}

template<typename TSeq>
inline void DataBase<TSeq>::add_user_data(
    std::vector< epiworld_double > x
)
{

    user_data.add(x);

}

template<typename TSeq>
inline void DataBase<TSeq>::add_user_data(
    epiworld_fast_uint k,
    epiworld_double x
)
{

    user_data.add(k, x);

}

template<typename TSeq>
inline UserData<TSeq> & DataBase<TSeq>::get_user_data()
{
    return user_data;
}

template<typename TSeq>
inline MapVec_type<int,int> DataBase<TSeq>::reproductive_number()
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

template<typename TSeq>
inline void DataBase<TSeq>::reproductive_number(
    std::string fn
) const {


    auto map = reproductive_number();

    std::ofstream fn_file(fn, std::ios_base::out);

    // Repeat the same error if the file doesn't exists
    if (!fn_file)
    {
        throw std::runtime_error(
            "Could not open file \"" + fn +
            "\" for writing.")
            ;
    }

    fn_file << 
        #ifdef EPI_DEBUG
        "thread " <<
        #endif
        "virus_id virus source source_exposure_date rt\n";


    for (auto & m : map)
        fn_file <<
            #ifdef EPI_DEBUG
            EPI_GET_THREAD_ID() << " " <<
            #endif
            m.first[0u] << " \"" <<
            virus_name[m.first[0u]] << "\" " <<
            m.first[1u] << " " <<
            m.first[2u] << " " <<
            m.second << "\n";

    return;

}

template<typename TSeq>
inline std::vector< epiworld_double > DataBase<TSeq>::transition_probability(
    bool print,
    bool normalize
) const {

    auto states_labels = model->get_states();
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

template<>
inline bool DataBase<std::vector<int>>::operator==(const DataBase<std::vector<int>> & other) const
{
    VECT_MATCH(
        virus_name, other.virus_name,
        "DataBase:: virus_name don't match"
        )

    EPI_DEBUG_FAIL_AT_TRUE(
        virus_sequence.size() != other.virus_sequence.size(),
        "DataBase:: virus_sequence don't match."
        )

    for (size_t i = 0u; i < virus_sequence.size(); ++i)
    {
        VECT_MATCH(
            virus_sequence[i], other.virus_sequence[i],
            "DataBase:: virus_sequence[i] don't match"
            )
    }

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
        tool_sequence,
        other.tool_sequence,
        "DataBase:: tool_sequence[i] don't match"
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
        "DataBase:: today_virus don't match."
        )
    
    for (size_t i = 0u; i < today_virus.size(); ++i)
    {
        VECT_MATCH(
            today_virus[i], other.today_virus[i],
            "DataBase:: today_virus[i] don't match"
            )
    }

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    if (today_tool.size() != other.today_tool.size())
        return false;
    
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
        "DataBase:: today_total don't match"
        )

    // Totals
    EPI_DEBUG_FAIL_AT_TRUE(
        today_total_nviruses_active != other.today_total_nviruses_active,
        "DataBase:: today_total_nviruses_active don't match."
        )

    // Transmission network
    VECT_MATCH(
        transmission_date,
        other.transmission_date,                 ///< Date of the transmission eve,
        "DataBase:: transmission_date[i] don't match"
        )

    VECT_MATCH(
        transmission_source,
        other.transmission_source,               ///< Id of the sour,
        "DataBase:: transmission_source[i] don't match"
        )

    VECT_MATCH(
        transmission_target,
        other.transmission_target,               ///< Id of the targ,
        "DataBase:: transmission_target[i] don't match"
        )

    VECT_MATCH(
        transmission_virus,
        other.transmission_virus,              ///< Id of the varia,
        "DataBase:: transmission_virus[i] don't match"
        )

    VECT_MATCH(
        transmission_source_exposure_date,
        other.transmission_source_exposure_date, ///< Date when the source acquired the varia,
        "DataBase:: transmission_source_exposure_date[i] don't match"
        )


    VECT_MATCH(
        transition_matrix,
        other.transition_matrix,
        "DataBase:: transition_matrix[i] don't match"
        )


    return true;

}

template<typename TSeq>
inline bool DataBase<TSeq>::operator==(const DataBase<TSeq> & other) const
{
    VECT_MATCH(
        virus_name,
        other.virus_name,
        "DataBase:: virus_name[i] don't match"
    )

    VECT_MATCH(
        virus_sequence,
        other.virus_sequence,
        "DataBase:: virus_sequence[i] don't match"
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
        tool_sequence,
        other.tool_sequence,
        "DataBase:: tool_sequence[i] don't match"
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
        transition_matrix,
        other.transition_matrix,
        "DataBase:: transition_matrix[i] don't match"
    )

    return true;

}

template<typename TSeq>
inline void DataBase<TSeq>::generation_time(
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

template<typename TSeq>
inline void DataBase<TSeq>::generation_time(
    std::string fn
) const
{

    std::vector< int > agent_id;
    std::vector< int > virus_id;
    std::vector< int > time;
    std::vector< int > gentime;

    generation_time(agent_id, virus_id, time, gentime);

    std::ofstream fn_file(fn, std::ios_base::out);

    // Throw an error if the file doesn't exists using throw
    if (!fn_file)
    {
        throw std::runtime_error(
            "DataBase::generation_time: "
            "Cannot open file " + fn + "."
        );
    }



    fn_file << 
        #ifdef EPI_DEBUG
        "thread " <<
        #endif
        "virus source source_exposure_date gentime\n";

    size_t n = agent_id.size();
    for (size_t i = 0u; i < n; ++i)
        fn_file <<
            #ifdef EPI_DEBUG
            EPI_GET_THREAD_ID() << " " <<
            #endif
            virus_id[i] << " " <<
            agent_id[i] << " " <<
            time[i] << " " <<
            gentime[i] << "\n";

    return;

}

#undef VECT_MATCH

#endif