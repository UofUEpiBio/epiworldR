#ifndef EPIWORLD_HOSPITALIZATIONSTRACKER_MEAT_HPP
#define EPIWORLD_HOSPITALIZATIONSTRACKER_MEAT_HPP


template<typename TSeq>
inline void HospitalizationsTracker<TSeq>::reset()
{
    _date.clear();
    _virus_id.clear();
    _tool_id.clear();
    _weight.clear();
}

template<typename TSeq>
inline void HospitalizationsTracker<TSeq>::record(
    Agent<TSeq> & agent,
    Model<TSeq> & model
)
{
    int current_date = model.today();
    
    // Get the virus ID (-1 if no virus)
    int v_id = -1;
    auto & virus = agent.get_virus();
    if (virus != nullptr)
        v_id = virus->get_id();
    
    // Get the number of tools
    size_t n_tools = agent.get_n_tools();
    
    if (n_tools == 0u)
    {
        // No tools: single record with tool_id = -1 and weight = 1.0
        _date.push_back(current_date);
        _virus_id.push_back(v_id);
        _tool_id.push_back(-1);
        _weight.push_back(1.0);
    }
    else
    {
        // Multiple tools: one record per tool with weight = 1/N
        double w = 1.0 / static_cast<double>(n_tools);
        for (size_t i = 0u; i < n_tools; ++i)
        {
            _date.push_back(current_date);
            _virus_id.push_back(v_id);
            _tool_id.push_back(agent.get_tool(static_cast<int>(i))->get_id());
            _weight.push_back(w);
        }
    }
}

template<typename TSeq>
inline void HospitalizationsTracker<TSeq>::get(
    int ndays,
    std::vector<int> & date,
    std::vector<int> & virus_id,
    std::vector<int> & tool_id,
    std::vector<int> & count,
    std::vector<double> & weight
) const
{
    // Clear output vectors
    date.clear();
    virus_id.clear();
    tool_id.clear();
    count.clear();
    weight.clear();
    
    if (ndays <= 0)
        return;
    
    // First, aggregate by (date, virus_id, tool_id) to get count and sum of weights
    // Key: (date, virus_id, tool_id), Value: (count, sum of weight)
    std::map<std::tuple<int, int, int>, std::pair<int, double>> aggregated;
    
    // Collect unique (virus_id, tool_id) combinations
    std::set<std::pair<int, int>> unique_combinations;
    
    for (size_t i = 0u; i < _date.size(); ++i)
    {
        auto key = std::make_tuple(_date[i], _virus_id[i], _tool_id[i]);
        aggregated[key].first += 1;  // Count
        aggregated[key].second += _weight[i];  // Weight
        unique_combinations.insert(std::make_pair(_virus_id[i], _tool_id[i]));
    }
    
    // If no hospitalizations occurred, still output the full time series 
    // with a single combination (virus_id=-1, tool_id=-1) and zeros
    if (unique_combinations.empty())
    {
        unique_combinations.insert(std::make_pair(-1, -1));
    }
    
    // For each unique (virus_id, tool_id) combination, output all days
    for (const auto & combo : unique_combinations)
    {
        int v_id = combo.first;
        int t_id = combo.second;
        
        for (int d = 0; d < ndays; ++d)
        {
            date.push_back(d);
            virus_id.push_back(v_id);
            tool_id.push_back(t_id);
            
            // Look up the count and weight for this (date, virus_id, tool_id) combination
            auto key = std::make_tuple(d, v_id, t_id);
            auto it = aggregated.find(key);
            if (it != aggregated.end())
            {
                count.push_back(it->second.first);
                weight.push_back(it->second.second);
            }
            else
            {
                count.push_back(0);
                weight.push_back(0.0);
            }
        }
    }
}

template<typename TSeq>
inline size_t HospitalizationsTracker<TSeq>::size() const
{
    return _date.size();
}

#endif
