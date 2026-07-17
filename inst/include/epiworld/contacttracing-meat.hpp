#ifndef EPIWORLD_CONTACTTRACING_MEAT_H
#define EPIWORLD_CONTACTTRACING_MEAT_H

#include "contacttracing-bones.hpp"

inline size_t ContactTracing::get_location(size_t row, size_t col) const
{
    return col * n_agents + row;
}

inline ContactTracing::ContactTracing()
{
    n_agents = 0u;
    max_contacts = 0u;
}

inline ContactTracing::ContactTracing(size_t n_agents, size_t max_contacts)
{
    this->n_agents = n_agents;
    this->max_contacts = max_contacts;

    contact_matrix.resize(n_agents * max_contacts, 0u);
    contacts_per_agent.resize(n_agents, 0);
    contact_date.resize(n_agents * max_contacts, 0);

    cached_contacts.resize(n_agents);
    up_to_date.assign(n_agents, false);
}

inline void ContactTracing::add_contact(size_t agent_a, size_t agent_b, size_t day)
{
    
    // Checking overflow
    size_t col_location = contacts_per_agent[agent_a] % max_contacts;
    size_t array_location = get_location(agent_a, col_location);

    contact_matrix[array_location] = agent_b;
    contact_date[array_location] = day;

    contacts_per_agent[agent_a] += 1;

    // Invalidate the cache for this agent
    up_to_date[agent_a] = false;

}

inline size_t ContactTracing::get_n_contacts(size_t agent)
{
    return contacts_per_agent[agent];
}

inline size_t ContactTracing::get_max_contacts() const
{
    return max_contacts;
}

inline std::pair< size_t, int> ContactTracing::get_contact(size_t agent, size_t idx)
{
    if (
        (idx >= contacts_per_agent[agent]) ||
        (idx >= max_contacts)
    )
        throw std::out_of_range("Index out of range in get_contact");

    size_t col_location = idx % max_contacts;
    size_t array_location = get_location(agent, col_location);

    return { contact_matrix[array_location], contact_date[array_location] };
}

inline void ContactTracing::reset(size_t n_agents, size_t max_contacts)
{

    this->n_agents = n_agents;
    this->max_contacts = max_contacts;

    contact_matrix.assign(n_agents * max_contacts, 0u);
    contacts_per_agent.assign(n_agents, 0u);
    contact_date.assign(n_agents * max_contacts, 0u);

    cached_contacts.assign(n_agents, std::vector<ContactRecord>());
    up_to_date.assign(n_agents, false);
}

inline const std::vector<ContactRecord> & ContactTracing::get_contacts(size_t agent)
{
    if (!up_to_date[agent])
    {
        // Rebuild the cache for this agent by grouping stored contacts by
        // contact id and collecting all recorded days into a std::set<int>.
        std::map<size_t, std::set<int>> contact_map;

        size_t actual_n = contacts_per_agent[agent];
        if (actual_n > max_contacts)
            actual_n = max_contacts;

        for (size_t i = 0u; i < actual_n; ++i)
        {
            size_t array_location = get_location(agent, i);
            contact_map[contact_matrix[array_location]].insert(
                static_cast<int>(contact_date[array_location])
            );
        }

        cached_contacts[agent].clear();
        cached_contacts[agent].reserve(contact_map.size());
        for (const auto & kv : contact_map)
            cached_contacts[agent].emplace_back(kv.first, kv.second);

        up_to_date[agent] = true;
    }

    return cached_contacts[agent];
}

inline void ContactTracing::print(size_t agent)
{

    size_t n_contacts = contacts_per_agent[agent];
    if (n_contacts > max_contacts)
        n_contacts = max_contacts;

    printf_epiworld("Agent %zu has %zu contacts: ", agent, n_contacts);
    for (size_t i = 0u; i < n_contacts; ++i)
    {
        auto [contact_id, contact_day] = get_contact(agent, i);
        printf_epiworld("(%zu, day %i) ", contact_id, contact_day);
    }
    printf_epiworld("\n");

    return;

}

#endif