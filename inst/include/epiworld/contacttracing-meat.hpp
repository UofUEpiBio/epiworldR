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
}

inline void ContactTracing::add_contact(size_t agent_a, size_t agent_b, size_t day)
{
    
    // Checking overflow
    size_t col_location = contacts_per_agent[agent_a] % max_contacts;
    size_t array_location = get_location(agent_a, col_location);

    contact_matrix[array_location] = agent_b;
    contact_date[array_location] = day;

    contacts_per_agent[agent_a] += 1;

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