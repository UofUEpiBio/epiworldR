#ifndef EPIWORLD_CONTACTLAYER_MEAT_HPP
#define EPIWORLD_CONTACTLAYER_MEAT_HPP

/**
 * @file contactlayer-meat.hpp
 * @brief Who takes part in a contact layer: its excluded statuses and its
 * closures (see `ContactLayer`).
 */

inline void ContactLayer::exclude_status(std::string status)
{

    if (in_run)
        throw std::logic_error(
            "The contact layer '" + name + "': the statuses it excludes are "
            "set between runs (close() and reopen() act during a run)."
        );

    if (
        std::find(excluded_statuses.begin(), excluded_statuses.end(), status) ==
        excluded_statuses.end()
    )
        excluded_statuses.push_back(std::move(status));

}

inline void ContactLayer::close()
{
    closed_next = true;
    if (!in_run)
        closed_base = closed = true;
}

inline void ContactLayer::reopen()
{
    closed_next = false;
    if (!in_run)
        closed_base = closed = false;
}

inline bool ContactLayer::excludes(epiworld_fast_uint state) const
{
    return (!model->state_excluded.empty() && model->state_excluded[state]) ||
        (!state_out.empty() && state_out[state]);
}

inline bool ContactLayer::excludes_agent(size_t agent_id) const
{
    return excludes(model->agent_state[agent_id]);
}

inline const std::vector< bool > * ContactLayer::silent() const
{
    if (!state_silent.empty())
        return &state_silent;

    return model->state_not_infectious.empty() ?
        nullptr : &model->state_not_infectious;
}

inline void ContactLayer::reset_participation(const Model & m)
{

    closed = closed_next = closed_base;
    closed_spans.clear();
    if (closed_base)
        closed_spans.push_back({0, INT_MAX});

    state_out.clear();
    state_silent.clear();
    if (excluded_statuses.empty())
        return;

    // The statuses this layer excludes, by id
    const std::vector< std::string > names = (m.get_statuses() != nullptr) ?
        m.get_statuses()->get_stages() : std::vector< std::string >{"Free"};

    std::vector< char > status_out(names.size(), 0);
    for (const auto & s : excluded_statuses)
    {
        auto it = std::find(names.begin(), names.end(), s);
        if (it == names.end())
        {
            std::string have;
            for (const auto & n : names)
                have += (have.empty() ? "\"" : ", \"") + n + "\"";

            throw std::invalid_argument(
                "The contact layer '" + name + "' excludes the status \"" + s +
                "\", which the model does not have (its statuses: " + have +
                ")."
            );
        }
        status_out[static_cast< size_t >(it - names.begin())] = 1;
    }

    // The engine states in those statuses
    const size_t ns = static_cast< size_t >(m.nstates);
    state_out.assign(ns, 0);
    for (size_t c = 0u; c < ns; ++c)
        if (status_out[static_cast< size_t >(m.nh_status_of(static_cast< unsigned int >(c)))])
            state_out[c] = 1;

    // What transmits nothing here: the model's (which includes the states
    // out of every layer), plus those out of this one
    state_silent = m.state_not_infectious;
    state_silent.resize(ns, false);
    for (size_t c = 0u; c < ns; ++c)
        state_silent[c] = state_silent[c] || (state_out[c] != 0);

}

inline void ContactLayer::apply_closure(int today)
{

    if (closed_next == closed)
        return;

    if (closed_next)
        closed_spans.push_back({today + 1, INT_MAX});
    else
        closed_spans.back().second = today;

    closed = closed_next;

}

template< typename F >
inline void ContactLayer::open_days(int a, int b, F f) const
{

    // The closed spans are in order and do not overlap
    int from = a;
    for (const auto & span : closed_spans)
    {
        if (span.second < from)
            continue;

        if (span.first > b)
            break;

        if (span.first > from)
            f(from, span.first - 1);

        if (span.second >= b)
            return;

        from = span.second + 1;
    }

    if (from <= b)
        f(from, b);

}

#endif
