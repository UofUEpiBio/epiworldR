#ifndef EPIWORLD_QUEUE_MEAT_HPP
#define EPIWORLD_QUEUE_MEAT_HPP

// The member functions of `Queue` (see queue-bones.hpp)

inline void Queue::shift(size_t id, epiworld_fast_int n)
{

    if (n == 0)
        return;

    epiworld_fast_int before = active[id];
    active[id] += n;

    // The agent enters or leaves the queue only when its count crosses zero.
    // Then the count of queued agents changes, and so does the agent's bit in
    // `bits`: word id / 64 (id >> 6), bit id % 64 (id & 63). Setting and
    // clearing that bit here is what keeps `for_each_nonzero()` -- the ordered
    // walk over queued agents that replaces scanning the whole population --
    // in step with the counts.
    if ((before == 0) && (active[id] != 0))
    {
        n_in_queue++;
        bits[id >> 6] |= (uint64_t(1) << (id & 63u));
    }
    else if ((before != 0) && (active[id] == 0))
    {
        n_in_queue--;
        bits[id >> 6] &= ~(uint64_t(1) << (id & 63u));
    }

}

inline void Queue::shift_neighbors(size_t id, epiworld_fast_int n)
{

    for (const auto * net : model->network_layers)
        for (size_t j : net->ties[id])
            shift(j, n);

    for (const auto * mix : model->group_layers)
    {
        const size_t g = mix->groups[id];
        if (g == layer::Mixing::NO_GROUP)
            continue;

        for (size_t j : mix->members[g])
            if (j != id)
                shift(j, n);
    }

}

inline void Queue::operator+=(Agent * p)
{

    const size_t id = static_cast< size_t >(p->id);
    everyone[id]++;
    shift(id, 1);
    shift_neighbors(id, 1);

}

inline void Queue::operator-=(Agent * p)
{

    const size_t id = static_cast< size_t >(p->id);
    everyone[id]--;
    shift(id, -1);
    shift_neighbors(id, -1);

}

inline bool Queue::tracks(size_t a, size_t b) const
{

    // The counters are sized by reset(), i.e. when a run starts. Editing the
    // network before that -- while the model is still being set up -- has no
    // queue to keep in step, and the counts are built from the finished network
    // anyway.
    return everyone.size() > (a > b ? a : b);

}

inline void Queue::notify_ties_changed(size_t a, int da, size_t b, int db)
{

    if (!tracks(a, b))
        return;

    shift(b, da * everyone[a]);
    shift(a, db * everyone[b]);

}

inline void Queue::notify_group_changed(
    size_t a,
    const std::vector< size_t > * from,
    const std::vector< size_t > * to
)
{

    if (!tracks(a, a))
        return;

    if (from != nullptr)
        for (size_t j : *from)
            if (j != a)
                notify_ties_changed(a, -1, j, -1);

    if (to != nullptr)
        for (size_t j : *to)
            if (j != a)
                notify_ties_changed(a, 1, j, 1);

}

inline void Queue::notify_edges_swapped(
    size_t a,
    size_t b,
    size_t c,
    size_t d,
    bool directed
)
{

    if (!tracks(a, b) || !tracks(c, d))
        return;

    // `a`'s entry for `b` now names `d`, and `c`'s entry for `d` now names `b`.
    epiworld_fast_int delta = everyone[c] - everyone[a];
    shift(b, delta);
    shift(d, -delta);

    if (directed)
        return;

    // `b`'s entry for `a` now names `c`, and `d`'s entry for `c` now names `a`.
    delta = everyone[d] - everyone[b];
    shift(a, delta);
    shift(c, -delta);

}

inline epiworld_fast_int Queue::operator[](epiworld_fast_uint i) const
{
    return active[i];
}

template< typename F >
inline void Queue::for_each_nonzero(F && f)
{

    const size_t nwords = bits.size();
    for (size_t w = 0u; w < nwords; ++w)
    {

        uint64_t word = bits[w];
        while (word != 0u)
        {

            unsigned int b = epi_ctz64(word);
            f((w << 6) + b);

            // Re-read the word: `f` may have queued or dequeued agents ahead
            // of this one. Only the bits above `b` are still to be visited.
            word = (b == 63u) ? 0u : (bits[w] & (~uint64_t(0) << (b + 1u)));

        }

    }

}

inline void Queue::reset()
{

    // Cleared unconditionally: counts that never went through `shift()` (or
    // that went negative) would otherwise survive into the next run.
    size_t n = model->size();
    active.assign(n, 0);
    everyone.assign(n, 0);
    bits.assign((n + 63u) / 64u, 0u);
    n_in_queue = 0;

}

inline bool Queue::operator==(const Queue & other) const 
{
    if (active.size() != other.active.size())
        return false;

    for (size_t i = 0u; i < active.size(); ++i)
    {
        if (active[i] != other.active[i])
            return false;
    }

    if (everyone.size() != other.everyone.size())
        return false;

    for (size_t i = 0u; i < everyone.size(); ++i)
    {
        if (everyone[i] != other.everyone[i])
            return false;
    }

    return true;
}

#endif
