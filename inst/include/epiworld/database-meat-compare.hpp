#ifndef EPIWORLD_DATABASE_MEAT_COMPARE_HPP
#define EPIWORLD_DATABASE_MEAT_COMPARE_HPP

/**
 * @file database-meat-compare.hpp
 * @brief Comparing databases (each part compares its own fields; with
 * `EPI_DEBUG`, a mismatch throws naming the field).
 */

#define VECT_MATCH(a, b, c) \
    EPI_DEBUG_FAIL_AT_TRUE(a.size() != b.size(), c) \
    for (size_t __i = 0u; __i < a.size(); ++__i) \
    {\
        EPI_DEBUG_FAIL_AT_TRUE(a[__i] != b[__i], c) \
    }

inline bool DataRegistry::operator==(const DataRegistry & other) const
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

    

    return true;

}

inline bool DiseaseRecord::operator==(const DiseaseRecord & other) const
{
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

inline bool DataBase::operator==(const DataBase & other) const
{

    if (!(reg == other.reg))
        return false;

    EPI_DEBUG_FAIL_AT_TRUE(
        sampling_freq != other.sampling_freq,
        "DataBase:: sampling_freq don't match."
    )

    return rec == other.rec;

}

#undef VECT_MATCH

#endif
