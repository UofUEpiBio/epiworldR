#ifndef EPIWORLD_CONTACTTRACING_BONES_H
#define EPIWORLD_CONTACTTRACING_BONES_H

#include <vector>
#include <stdexcept>
#include "config.hpp"

/** 
 * @brief Class for tracing contacts between agents
 * @details
 * The class assumes that contacts are stored in a matrix-like
 * structure, where rows are agents and columns are contacts in a 
 * column-major format. Each entry stores the id of the contact
 * and the day of the contact is stored in a parallel array.
 * 
 * The maximum number of contacts per agent is fixed at initialization
 * time. If more contacts are added, they will overwrite previous
 * contacts in a circular buffer fashion.
 * */
class ContactTracing
{
private:

    std::vector< size_t > contact_matrix;
    std::vector< size_t > contacts_per_agent;
    std::vector< size_t > contact_date;

    size_t n_agents;
    size_t max_contacts;

    size_t get_location(size_t row, size_t col) const;

public:

    ContactTracing();

    /**
     * @brief Construct a new Contact Tracing object
     * 
     * @param n_agents Agents in the system, usually `Model<TSeq>::size()`.
     * @param max_contacts Maximum number of contacts to track per agent.
     */
    ContactTracing(size_t n_agents, size_t max_contacts);

    /**
     * @brief Add a contact between two agents at a given day.
     * 
     * @param agent_a Agent id (usually infectious agent)
     * @param agent_b Agent id (usually susceptible agent)
     * @param day Day of the contact (usually Model<TSeq>::today()).
     */
    void add_contact(size_t agent_a, size_t agent_b, size_t day);

    /**
     * @brief Get the number of contacts for an agent
     * 
     * @param agent Agent id
     * @return size_t Number of contacts recorded for that agent (can be more than max_contacts)
     */
    size_t get_n_contacts(size_t agent); 
    
    /**
     * @brief Get the contact object
     * 
     * @param agent Source agent (id)
     * @param idx Index of the contact (0 to get_n_contacts(agent)-1)
     * @return std::pair<size_t, size_t> with (contact_id, contact_day)
     * @throws std::out_of_range if idx is out of range.
     */
    std::pair<size_t, size_t> get_contact(size_t agent, size_t idx);

    /**
     * @brief Reset the contact tracing data
     * 
     * Usually called by `Model<TSeq>::reset()`.
     * 
     * @param n_agents Number of agents
     * @param max_contacts Maximum number of contacts to track per agent
     */
    void reset(
        size_t n_agents,
        size_t max_contacts
    );

    /**
     * @brief Print the contacts of an agent 
     * 
     * @param agent Agent id
     */
    void print(size_t agent);
};

#endif
