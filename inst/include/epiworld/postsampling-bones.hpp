#ifndef EPIWORLD_POSTSAMPLING_BONES_HPP
#define EPIWORLD_POSTSAMPLING_BONES_HPP

#include <cstdint>
#include <cstddef>
#include <functional>
#include <vector>

template<typename TSeq> class Agent;
template<typename TSeq> class Model;

/**
 * @brief Read-only view of the agents one infectious agent was in contact with.
 *
 * @details Passed to a `PostSamplingFun`. It points into the model's scratch
 * storage, so it is valid only during the callback: do not keep it (copy the
 * ids if you need them later). The contacts are a multiset: repeated contacts
 * appear repeatedly, and contacts that did not transmit are included. The
 * order of the contacts within the view is unspecified.
 */
class SampledContactsView {
private:
    const uint32_t * first = nullptr;
    size_t n = 0u;
public:
    SampledContactsView() = default;
    SampledContactsView(const uint32_t * first_, size_t n_) : first(first_), n(n_) {}
    const uint32_t * begin() const { return first; }
    const uint32_t * end() const { return first + n; }
    size_t size() const { return n; }
    bool empty() const { return n == 0u; }
    size_t operator[](size_t i) const { return first[i]; }
};

/**
 * @brief Callback run after the contacts of a step were sampled.
 *
 * @details Called once per infectious agent that had at least one sampled
 * contact in the step, in ascending order of the agent's id, after all the
 * susceptible agents were updated and right before the events are applied.
 * The agents still have the state they had when the contacts were sampled;
 * events the callback queues are applied in the same step.
 *
 * Arguments: the infectious agent, the ids of the agents it was in contact
 * with (see `SampledContactsView`), and the model.
 */
template<typename TSeq>
using PostSamplingFun = std::function<
    void(Agent<TSeq> *, const SampledContactsView &, Model<TSeq> *)
>;

/**
 * @brief Model-local storage used to batch sampled contacts.
 *
 * @details Pairs `(infectious, contacted)` are appended while sampling and
 * grouped by infectious agent once, at the end of the sampling phase. Lengths
 * are cleared, not released, between steps. Copies of a model start with empty
 * scratch (nothing is shared or copied); the model sizes `counts` for its
 * population the first time it needs it (see `Model::post_sampling_dispatch()`).
 */
struct PostSamplingScratch {

    std::vector< uint32_t > pairs;    ///< Flat: infectious, contacted, ...
    std::vector< uint32_t > grouped;  ///< Contacted ids, by infectious agent
    std::vector< uint32_t > touched;  ///< Distinct infectious ids (sorted)
    std::vector< size_t > starts;     ///< [touched] Where the batch starts
    std::vector< uint32_t > counts;   ///< [agent] Zero between steps

    PostSamplingScratch() = default;
    PostSamplingScratch(const PostSamplingScratch &) {}
    PostSamplingScratch & operator=(const PostSamplingScratch &) { reset(0u); return *this; }

    /// Empties the batch and restores the invariant of `counts` (all zero).
    void clear()
    {
        for (auto id : touched)
            counts[id] = 0u;
        pairs.clear();
        grouped.clear();
        touched.clear();
        starts.clear();
    }

    /// Sizes the storage for a population (call at reset).
    void reset(size_t n_agents)
    {
        pairs.clear();
        grouped.clear();
        touched.clear();
        starts.clear();
        counts.assign(n_agents, 0u);
    }

};

#endif
