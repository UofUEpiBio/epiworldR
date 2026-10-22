#ifndef EPIWORLD_QUEUE_BONES_HPP
#define EPIWORLD_QUEUE_BONES_HPP

/**
 * @brief Controls which agents are verified at each step
 * 
 * @details The idea is that only agents who are either in
 * an infected state or have an infected neighbor should be
 * checked. Otherwise it makes no sense (no chance to recover
 * or capture the disease).
 * 
 * **Implementation details:**
 * <a href="../impl/queueing-system.md">Queueing System</a>
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Queue
{
    friend class Model<TSeq>;

private:

    /**
     * @brief Count of ego's neighbors in queue (including ego)
     */
    std::vector< epiworld_fast_int > active;

    /**
     * @brief Outstanding `Everyone` registrations per agent.
     *
     * @details `active[i]` is `everyone[i]` plus the sum of `everyone[j]` over
     * the agents `j` whose neighbor lists contain `i`, in every network layer
     * of the model (`operator+=` credits the agents in the registered agent's
     * own lists). In an undirected network those are `i`'s neighbors; in a
     * directed one they are the agents that list `i`. So this is what a tie is worth to the far end of
     * it: adding an edge next to a registered agent is a `+everyone[]`
     * over there, and removing one is the mirror (see `notify_edge_added`).
     *
     * Keeping the count explicit -- rather than inferring "is this agent
     * contributing?" from whether it carries a virus -- means the bookkeeping
     * stays right no matter what registered the agent: a virus, a tool, or an
     * explicit queue argument to `Agent::change_state`.
     */
    std::vector< epiworld_fast_int > everyone;

    /**
     * @brief One bit per agent: bit `i` is set iff `active[i] != 0`.
     *
     * @details This is what lets `Model::update_state()` visit only the agents
     * in the queue, in ascending id order, without scanning the whole
     * population: it walks the set bits (see `for_each_nonzero`). Every change
     * to `active` goes through `shift()`, which keeps the two in step.
     */
    std::vector< uint64_t > bits;

    Model<TSeq> * model = nullptr;
    int n_in_queue = 0;

    /// @brief Adds `n` to `active[id]`, keeping `n_in_queue` and `bits` in step.
    void shift(size_t id, epiworld_fast_int n);

    /// @brief Whether the counters have been sized to cover these two agents.
    bool tracks(size_t a, size_t b) const;

    /// @brief Adds `n` to the counts of `id`'s neighbors in every network.
    void shift_neighbors(size_t id, epiworld_fast_int n);

    // Auxiliary variable that checks how many steps
    // left are there
    // int n_steps_left;
    // bool queuing_started   = false;

public:

    void operator+=(Agent<TSeq> * p);
    void operator-=(Agent<TSeq> * p);

    /**
     * @brief How many registrations cover agent `i` (itself or a neighbor).
     *
     * @details Read-only: the count must only change through `+=`, `-=` and
     * the edge notifications, which keep the ordered set of queued agents in
     * step with it.
     */
    epiworld_fast_int operator[](epiworld_fast_uint i) const;

    /**
     * @brief Calls `f(i)` for every agent with a non-zero count, in ascending
     * id order.
     *
     * @details Costs O(N / 64 + number of queued agents), instead of a scan of
     * the whole population. The count of an agent is read when the walk
     * reaches it, so if `f` changes the queue (e.g., a state function that
     * adds a tie), agents further ahead see the change -- exactly what a plain
     * `for (i = 0; i < N; ++i) if (queue[i] ...)` loop would do.
     */
    template< typename F >
    void for_each_nonzero(F && f);

    /**
     * @name Keep the queue in step with a change to the contact network
     *
     * @details The queue counts, for every agent, how many of its neighbors are
     * registered as active. That count is built when an agent is registered
     * (`operator+=`) and unwound when it is deregistered (`operator-=`), both
     * walking the agent's neighbors *as they are at that moment*. Changing the
     * network in between would leave the two walks disagreeing, and an agent
     * whose count drifted to zero is silently skipped by
     * `Model::update_state()`.
     *
     * This keeps the counts exact as the change happens, in constant time: an
     * agent that gains a tie to another hands it whatever it contributes, and
     * one that loses the tie takes it back.
     *
     * @param a,b The two ends of the tie that was just added or removed.
     * @param da,db How `a`'s list changed (`+1` gained `b`, `-1` lost it, `0`
     * unchanged), and `b`'s.
     */
    void notify_ties_changed(size_t a, int da, size_t b, int db);

    /**
     * @brief Keep the queue in step with `layer::Network::swap_ties()`.
     *
     * @details `a`'s tie to `b` and `c`'s tie to `d` traded ends: (a-b), (c-d)
     * became (a-d), (c-b). Nobody's degree changed, only whose registrations
     * reach whom: `b` now gets what `c` contributes instead of what `a` does,
     * and `d` the reverse. In an undirected network the other ends traded too,
     * so `a` gets what `d` contributes instead of `b`, and `c` the reverse.
     * Four shifts, whatever the degrees.
     *
     * In a directed network only `a`'s and `c`'s lists change, and the queue
     * follows each agent's own list (`operator+=`), so only `b`'s and `d`'s
     * counts change.
     *
     * @param a,c The agents whose neighbors were swapped.
     * @param b,d Their neighbors before the swap (`b` was `a`'s, `d` was `c`'s).
     * @param directed Whether only `a`'s and `c`'s lists changed.
     */
    void notify_edges_swapped(
        size_t a,
        size_t b,
        size_t c,
        size_t d,
        bool directed
    );

    // void initialize(Model<TSeq> * m, Agent<TSeq> * p);
    void reset();

    bool operator==(const Queue<TSeq> & other) const;
    bool operator!=(const Queue<TSeq> & other) const {return !operator==(other);};

    static const int NoOne    = 0;
    static const int OnlySelf = 1;
    static const int Everyone = 2;

};

template<typename TSeq>
inline void Queue<TSeq>::shift(size_t id, epiworld_fast_int n)
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

template<typename TSeq>
inline void Queue<TSeq>::shift_neighbors(size_t id, epiworld_fast_int n)
{

    for (const auto * net : model->network_layers)
        for (size_t j : net->ties[id])
            shift(j, n);

}

template<typename TSeq>
inline void Queue<TSeq>::operator+=(Agent<TSeq> * p)
{

    const size_t id = static_cast< size_t >(p->id);
    everyone[id]++;
    shift(id, 1);
    shift_neighbors(id, 1);

}

template<typename TSeq>
inline void Queue<TSeq>::operator-=(Agent<TSeq> * p)
{

    const size_t id = static_cast< size_t >(p->id);
    everyone[id]--;
    shift(id, -1);
    shift_neighbors(id, -1);

}

template<typename TSeq>
inline bool Queue<TSeq>::tracks(size_t a, size_t b) const
{

    // The counters are sized by reset(), i.e. when a run starts. Editing the
    // network before that -- while the model is still being set up -- has no
    // queue to keep in step, and the counts are built from the finished network
    // anyway.
    return everyone.size() > (a > b ? a : b);

}

template<typename TSeq>
inline void Queue<TSeq>::notify_ties_changed(size_t a, int da, size_t b, int db)
{

    if (!tracks(a, b))
        return;

    shift(b, da * everyone[a]);
    shift(a, db * everyone[b]);

}

template<typename TSeq>
inline void Queue<TSeq>::notify_edges_swapped(
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

template<typename TSeq>
inline epiworld_fast_int Queue<TSeq>::operator[](epiworld_fast_uint i) const
{
    return active[i];
}

template<typename TSeq>
template< typename F >
inline void Queue<TSeq>::for_each_nonzero(F && f)
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

template<typename TSeq>
inline void Queue<TSeq>::reset()
{

    // Cleared unconditionally: counts that never went through `shift()` (or
    // that went negative) would otherwise survive into the next run.
    size_t n = model->size();
    active.assign(n, 0);
    everyone.assign(n, 0);
    bits.assign((n + 63u) / 64u, 0u);
    n_in_queue = 0;

}

template<typename TSeq>
inline bool Queue<TSeq>::operator==(const Queue<TSeq> & other) const 
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