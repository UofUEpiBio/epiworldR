#ifndef EPIWORLD_CONTACTTRACING_BONES_H
#define EPIWORLD_CONTACTTRACING_BONES_H

#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include "config.hpp"

/**
 * @brief Represents a single contact relationship for an agent.
 * @details
 * Each instance holds the id of the contacted agent and the set of
 * simulation days on which the contact occurred.  Instances are
 * returned by `ContactTracing::get_contacts()`.
 */
class ContactRecord
{
private:
    size_t m_contact_id;
    std::set<int> m_times;

public:
    ContactRecord(size_t contact_id, std::set<int> times)
        : m_contact_id(contact_id), m_times(std::move(times)) {}

    /**
     * @brief Return the id of the contacted agent.
     */
    size_t get_contact_id() const { return m_contact_id; }

    /**
     * @brief Return the set of days on which the contact was recorded.
     * @return A const reference to the set of simulation days.
     */
    const std::set<int> & get_times() const { return m_times; }
};

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

    // Cache: for each agent, the list of unique contacts with their dates.
    // up_to_date[a] is false whenever add_contact(a, ...) is called, and
    // set to true after get_contacts(a) rebuilds the cache.
    std::vector< std::vector<ContactRecord> > cached_contacts;
    std::vector< bool > up_to_date;

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

    size_t get_max_contacts() const;
    
    /**
     * @brief Get the contact object
     * 
     * @param agent Source agent (id)
     * @param idx Index of the contact (0 to get_n_contacts(agent)-1)
     * @return std::pair<size_t, int> with (contact_id, contact_day)
     * @throws std::out_of_range if idx is out of range.
     */
    std::pair<size_t, int> get_contact(size_t agent, size_t idx);

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
     * @brief Get all unique contacts of an agent with the days they occurred.
     * @details
     * Returns a vector of `ContactRecord` objects, one per unique contacted
     * agent.  Each record exposes:
     *
     * - `get_contact_id()` – the id of the other agent.
     * - `get_times()` – the set of simulation days on which the contact was
     *   recorded in the current (non-reset) window.
     *
     * The result is lazily computed on the first call after any new contact is
     * recorded for this agent and cached until the next `add_contact()` call
     * for the same agent.
     *
     * @param agent Source agent id.
     * @return Const reference to the cached vector of ContactRecord objects.
     */
    const std::vector<ContactRecord> & get_contacts(size_t agent);

    /**
     * @brief Print the contacts of an agent 
     * 
     * @param agent Agent id
     */
    void print(size_t agent);
};

#endif
