#ifndef EPIWORLD_HOSPITALIZATIONSTRACKER_BONES_HPP
#define EPIWORLD_HOSPITALIZATIONSTRACKER_BONES_HPP

template<typename TSeq>
class Agent;

template<typename TSeq>
class Model;

/**
 * @brief Class to track hospitalizations in an epidemiological model.
 * 
 * @details
 * This class keeps track of hospitalizations in a model. For each hospitalized
 * agent, it records the date, virus ID, and tool IDs with appropriate weights.
 * 
 * Since agents always have at most one virus but may have multiple tools,
 * if an agent has N tools, then N records are created, each with weight
 * equal to 1/N. If the agent has no tools, a single record is created with
 * tool_id = -1 and weight = 1.0.
 * 
 * The `get()` method returns a summary of the hospitalizations grouped by
 * date, virus_id, and tool_id, with weight summed across all matching
 * records.
 * 
 * @tparam TSeq Type of sequence (should match `TSeq` across the model)
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class HospitalizationsTracker {

private:

    std::vector<int> _date;           ///< Date of hospitalization
    std::vector<int> _virus_id;       ///< ID of the virus causing hospitalization
    std::vector<int> _tool_id;        ///< ID of the tool (-1 if no tools)
    std::vector<double> _weight;      ///< Weight for the tool (1/N where N is tool count)

public:

    HospitalizationsTracker() = default;

    /**
     * @brief Reset the tracker by clearing all data.
     */
    void reset();

    /**
     * @brief Record a hospitalization event for an agent.
     * 
     * @param agent Reference to the agent being hospitalized.
     * @param model Reference to the model.
     * 
     * @details
     * For each hospitalization, the method records:
     * - The current date from the model
     * - The virus ID from the agent's virus
     * - For each tool the agent has, a separate record with weight = 1/N
     *   where N is the number of tools
     * - If the agent has no tools, a single record with tool_id = -1 and
     *   weight = 1.0
     */
    void record(Agent<TSeq> & agent, Model<TSeq> & model);

    /**
     * @brief Get the full time series of hospitalizations.
     * 
     * @param ndays Number of days in the simulation (0 to ndays-1).
     * @param date Output vector for dates.
     * @param virus_id Output vector for virus IDs.
     * @param tool_id Output vector for tool IDs.
     * @param count Output vector for counts (number of hospitalized individuals).
     * @param weight Output vector for summed weights (fractional contribution 
     *   based on tool distribution).
     * 
     * @details
     * Returns the full time series of hospitalization data. For each unique 
     * (virus_id, tool_id) combination observed, returns an entry for every 
     * day from 0 to ndays-1.
     * 
     * The `count` vector contains the actual number of individuals hospitalized 
     * for that (date, virus_id, tool_id) combination. This is useful for 
     * answering questions like "how many total people were hospitalized?" 
     * regardless of their tools.
     * 
     * The `weight` vector contains fractional contributions: if an agent has N 
     * tools, each tool gets weight = 1/N. Summing weights across all tool_ids 
     * for a given date and virus_id gives the total number of hospitalizations.
     */
    void get(
        int ndays,
        std::vector<int> & date,
        std::vector<int> & virus_id,
        std::vector<int> & tool_id,
        std::vector<int> & count,
        std::vector<double> & weight
    ) const;

    /**
     * @brief Get the number of raw records in the tracker.
     * @return Number of records.
     */
    size_t size() const;

};

#endif
