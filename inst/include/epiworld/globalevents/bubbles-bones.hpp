#ifndef EPIWORLD_GLOBALEVENTS_BUBBLES_BONES_HPP
#define EPIWORLD_GLOBALEVENTS_BUBBLES_BONES_HPP

// Standard library headers (vector, memory, string, unordered_map, utility,
// algorithm, stdexcept) are included at global scope by epiworld.hpp; this file
// is only ever included from within `namespace epiworld`, so it must not
// re-include system headers here.
#include "../config.hpp"

/**
 * @brief Flavor of the social-bubble intervention (how bubbles are formed).
 * @ingroup globalevents
 *
 * @details Both flavors produce a partition of the population into disjoint
 * bubbles in which households are never split, and both form bubbles along
 * *existing* contacts (see Bubbles).
 */
enum class BubbleFlavor {
    /**
     * Household-level rule: whole households bubble together, up to
     * `group_size` households per bubble. A bubble is grown from a seed
     * household by repeatedly absorbing a random household that is connected
     * to it in the contact network. `group_size == 1` leaves every household
     * on its own, i.e. a strict household-only lockdown.
     *
     * Models e.g. the Belgian rule of 10 May 2020 ("a household may form one
     * fixed bubble with other households").
     */
    Household,
    /**
     * Individual-level rule: each agent picks up to `group_size` peers among
     * its existing contacts outside its own household. Every pick merges the
     * two agents' households, and bubbles are the resulting connected groups
     * of households -- so a teenager choosing a partner brings both households
     * into one bubble.
     *
     * Models e.g. the Belgian rule of 19 October 2020 (`group_size == 1`, "one
     * close contact per person").
     */
    Peer
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
 * @details A "social bubble" policy lets people keep seeing a small, fixed set
 * of others while cutting off the rest of their contacts. `Bubbles` implements
 * this for models built on an explicit contact network (e.g. `ModelSEIR`,
 * `ModelSIR` after `agents_smallworld()` / `agents_from_edgelist()`), and is
 * designed for the household-based rules used during COVID-19.
 *
 * ## How it works
 *
 * The contact network is **not modified**. Instead, `deploy()` attaches a
 * `Tool` ("Social bubble") to every agent. When a susceptible agent `p` is
 * exposed to an infectious neighbor, the tool's susceptibility-reduction
 * function identifies the transmitter through the virus (`v->get_agent()`) and
 * compares the two agents' bubble labels:
 *
 * - different bubbles -> reduction `1.0`, transmission is blocked;
 * - same bubble       -> reduction `1 - within_factor`, i.e. transmission is
 *                        scaled by `within_factor` (social distancing);
 * - outside the policy window -> reduction `0.0`, no effect.
 *
 * Since reductions combine as `1 - prod(1 - r_i)`, a reduction of `1.0` zeroes
 * the transmission probability regardless of any other tools the agent carries.
 * Contact weights are uniform in these models, so suppressing transmission on
 * out-of-bubble contacts is equivalent to deleting those contacts, while
 * keeping the network intact for other purposes (contact tracing, output).
 *
 * ## Forming bubbles (why ties matter)
 *
 * Households are declared with a per-agent `household_id` vector (one entry per
 * agent, indexed by agent id). Bubbles are always:
 *
 * - **disjoint** -- every agent belongs to exactly one bubble;
 * - **household-preserving** -- a household is never split; and
 * - **connection-aware** -- bubbles only ever join households that are actually
 *   connected in the contact network.
 *
 * The last point is essential. Because the intervention can only suppress
 * transmission along existing edges and never creates new ones, putting two
 * households that share no contact into the same bubble changes nothing at all.
 * Pairing households at random would therefore leave `group_size` inert --
 * behaving like a strict lockdown no matter how large the bubbles are. It also
 * mirrors the real policy: a household chooses a bubble partner it already
 * socialises with. See `BubbleFlavor` for the two rules.
 *
 * A household with no unassigned connected partner simply ends up in a smaller
 * bubble (possibly alone), so the number of bubbles is at least
 * `ceil(n_households / group_size)`.
 *
 * ## Scheduling
 *
 * `start_day` is the first day the policy applies; `end_day` is the (exclusive)
 * day it lifts, or `< 0` to never lift. With `rewire_every > 0` the partition is
 * re-randomised every that-many days, modelling policies whose bubbles change
 * over time (e.g. contacts renewed weekly); `0` keeps a fixed bubble.
 *
 * The initial partition is computed at reset time via the tool's distribution
 * function, so it is recomputed for each replicate of `run_multiple()` using
 * that replicate's seed and is already in force on day 1. Re-randomisations are
 * applied by a daily global event and take effect the following step.
 *
 * ## Example
 *
 * ```cpp
 * epimodels::ModelSEIR<> model("flu", 0.01, 0.1, 4.5, 1.0/7.0);
 * model.agents_smallworld(10000, 8, false, 0.05);
 *
 * std::vector< size_t > household_id(10000);
 * for (size_t i = 0u; i < 10000; ++i)
 *     household_id[i] = i / 3;             // households of three
 *
 * // Two households per bubble from day 10, halving within-bubble transmission.
 * Bubbles<> bubbles(
 *     household_id, BubbleFlavor::Household,
 *     2,      // group_size
 *     0.5,    // within_factor
 *     10      // start_day
 * );
 * bubbles.deploy(model);
 *
 * model.run(100, 1231);
 * ```
 *
 * @note All copies of a `Bubbles` object share one `BubbleState`, so replicates
 * in `run_multiple()` must run on a single thread (`nthreads = 1`).
 *
 * @note Not yet modelled: fixed out-of-bubble "sport partners", travel between
 * regions, and caps on the number of individuals (as opposed to households) a
 * bubble may contain.
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
     * @brief Configure a social-bubble policy.
     *
     * @param household_id Household label of each agent, indexed by agent id.
     *        Its length must equal the number of agents in the model. Labels are
     *        arbitrary (they need not be consecutive); agents sharing a label
     *        form a household and are always placed in the same bubble.
     * @param flavor Whether bubbles are chosen by households or by individuals
     *        (see `BubbleFlavor`).
     * @param group_size For `Household`, the maximum number of households per
     *        bubble (`1` = strict household-only lockdown). For `Peer`, the
     *        maximum number of external peers each agent may pick.
     * @param within_factor Multiplier applied to transmission between agents of
     *        the same bubble, in `[0, 1]`: `1.0` leaves within-bubble
     *        transmission untouched, `0.5` halves it (physical distancing), and
     *        `0.0` blocks it entirely. Transmission between bubbles is always
     *        blocked.
     * @param start_day First day on which the policy applies.
     * @param end_day Day on which the policy is lifted (exclusive). Use a
     *        negative value for a policy that never ends. Must be greater than
     *        `start_day`.
     * @param rewire_every Re-randomise the bubbles every this-many days, for
     *        policies whose contacts change over time. `0` keeps the bubbles
     *        fixed for the whole intervention.
     * @param name Name given to the tool and to the scheduler event; useful to
     *        look them up on the model afterwards.
     *
     * @throws std::range_error if `within_factor` is outside `[0, 1]`, if
     *         `group_size` is zero for the `Household` flavor, or if `end_day`
     *         is non-negative and not greater than `start_day`.
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
     * @brief Install the intervention on a model.
     *
     * Adds the bubble `Tool` (distributed to every agent, which also computes
     * the bubble partition at each reset) and, when `rewire_every > 0`, the
     * global event that re-randomises the bubbles.
     *
     * Call this **after** the model's agents and contact network exist (e.g.
     * after `agents_smallworld()`), since the household grouping is derived from
     * the network, and before `run()`.
     *
     * @param model Model to install the intervention on.
     * @throws std::length_error if `household_id` does not have exactly one
     *         entry per agent.
     */
    void deploy(Model<TSeq> & model);

    /**
     * @brief Shared state of the intervention (bubble labels and bookkeeping).
     *
     * The state is shared by every copy of this object and by the model's tool
     * and event, and is refreshed on each run.
     */
    std::shared_ptr< BubbleState > get_state() const;

    /**
     * @brief Current bubble label of every agent, indexed by agent id.
     *
     * Two agents may transmit only when their labels are equal. Populated on the
     * first reset (i.e. once the model has been run); empty before then.
     */
    const std::vector< int > & get_bubble_id() const;

    /// @brief The rule used to form bubbles.
    BubbleFlavor get_flavor() const;

};

#endif
