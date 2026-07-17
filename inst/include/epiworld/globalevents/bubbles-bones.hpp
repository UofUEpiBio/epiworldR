#ifndef EPIWORLD_GLOBALEVENTS_BUBBLES_BONES_HPP
#define EPIWORLD_GLOBALEVENTS_BUBBLES_BONES_HPP

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include "../config.hpp"

/**
 * @brief Flavor of the social-bubble intervention.
 * @ingroup globalevents
 */
enum class BubbleFlavor {
    Household, ///< Households are grouped into bubbles of `group_size` households.
    Peer       ///< Each agent picks up to `group_size` external peers; picks merge households.
};

/**
 * @brief Shared, mutable state of a Bubbles intervention.
 *
 * Held via `std::shared_ptr` and captured by both the bubble `Tool`
 * (read-only) and the scheduler `GlobalEvent` (read-write), so the two stay in
 * sync for the lifetime of the model. The bubble partition is a per-agent
 * integer label (`bubble_id`): two agents transmit only when their labels match.
 */
struct BubbleState {
    std::vector< int > bubble_id;  ///< Per-agent bubble label (index = agent id); -1 = unassigned.
    int last_sim_id = -1;          ///< Sim id the current partition was computed for.
    int last_epoch  = -1;          ///< Rewiring epoch the current partition was computed for.
};

/**
 * @brief Social-bubble contact-restriction intervention for network models.
 * @ingroup globalevents
 *
 * @details Confines each agent's *effective* contacts to a small, disjoint
 * "bubble" during an outbreak, without editing the contact network. It attaches
 * a `Tool` ("Social bubble") to every agent whose susceptibility-reduction
 * function inspects the infectious neighbor (`v->get_agent()`): if the two
 * agents are in different bubbles the reduction is `1.0` (transmission blocked);
 * if they share a bubble it is `1 - within_factor` (social distancing). Because
 * contact weights are uniform in the network models, suppressing transmission on
 * out-of-bubble contacts is equivalent to removing those contacts.
 *
 * Households are identified by a per-agent `household_id` vector (one entry per
 * agent). Two flavors form the partition:
 * - `Household`: households are shuffled and chunked into groups of `group_size`
 *   households each (`group_size == 1` is the household-only lockdown baseline).
 * - `Peer`: each agent picks up to `group_size` peers among its existing external
 *   neighbors (a neighbor in a different household); each pick merges the two
 *   households, and bubbles are the resulting connected components.
 *
 * Scheduling is driven by `start_day`, `end_day` (`< 0` = never ends) and
 * `rewire_every` (`> 0` re-randomizes the partition every that-many days, for
 * "changing" bubble policies). The partition is (re)computed by a daily scheduler
 * event and takes effect the following simulation step.
 *
 * @note Uses a single shared `BubbleState` per model, so parallel replicates in
 * `run_multiple` must use a single thread.
 *
 * @tparam TSeq Sequence type (should match `TSeq` across the model).
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class Bubbles {
private:

    std::vector< size_t > household_id;
    BubbleFlavor flavor;
    size_t group_size;               ///< households per bubble (Household) or max peers (Peer).
    epiworld_double within_factor;   ///< within-bubble transmission multiplier in [0, 1].
    int start_day;
    int end_day;
    int rewire_every;
    std::string name;
    std::shared_ptr< BubbleState > state;

    void compute_partition(Model<TSeq> * model) const;
    void partition_household(Model<TSeq> * model) const;
    void partition_peer(Model<TSeq> * model) const;

public:

    /**
     * @param household_id Per-agent household label (size must equal the number
     *        of agents).
     * @param flavor `BubbleFlavor::Household` or `BubbleFlavor::Peer`.
     * @param group_size Households per bubble (Household) or maximum external
     *        peers per agent (Peer).
     * @param within_factor Within-bubble transmission multiplier (1.0 = no
     *        distancing, 0.0 = fully blocked). Must be in [0, 1].
     * @param start_day First day the policy is in effect.
     * @param end_day Day the policy ends (exclusive); `< 0` means never.
     * @param rewire_every Re-randomize the partition every this-many days;
     *        `<= 0` keeps a fixed partition.
     * @param name Name used for the tool and scheduler event.
     */
    Bubbles(
        std::vector< size_t > household_id,
        BubbleFlavor flavor,
        size_t group_size,
        epiworld_double within_factor = 1.0,
        int start_day = 0,
        int end_day = -1,
        int rewire_every = 0,
        std::string name = "Social bubble"
    );

    /**
     * @brief Install the bubble tool and scheduler event on a model.
     * Must be called after the model's agents/network are built.
     */
    void deploy(Model<TSeq> & model);

    std::shared_ptr< BubbleState > get_state() const;
    const std::vector< int > & get_bubble_id() const;
    BubbleFlavor get_flavor() const;

};

#endif
