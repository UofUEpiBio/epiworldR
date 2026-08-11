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
     * Individual-level rule: each agent nominates up to `group_size` peers
     * among its existing contacts outside its own household. A nomination
     * merges the two agents' households -- so a teenager choosing a partner
     * brings both households into one bubble.
     *
     * Nominations are only accepted while the bubble stays within
     * `max_households` households; this cap is what enforces the policy's
     * exclusivity. Without it the merges percolate: with several members per
     * household each nominating someone, the household graph becomes connected
     * and everyone ends up in a single bubble, imposing no restriction at all.
     *
     * Models e.g. the Belgian rule of 19 October 2020 (`group_size == 1`, "one
     * close contact per person").
     */
    Peer
};

template<typename TSeq>
class Bubbles;

/**
 * @brief Tool carried by every agent under a Bubbles policy.
 * @ingroup globalevents
 *
 * @details The tool does not hold the bubble partition -- it holds a pointer to
 * the `Bubbles` intervention that owns it, resolved by name from the model the
 * first time the tool is used. `clone_ptr()` clears that pointer, exactly as
 * `ToolVaccine` clears its per-agent immunity, so a copy of the tool always
 * re-resolves against the model it ends up in. Copies of a model therefore
 * never reach into the model they were copied from, which is what lets bubble
 * models run in parallel replicates.
 *
 * Users do not instantiate this directly; `Bubbles::deploy()` does.
 *
 * @tparam TSeq Sequence type (should match `TSeq` across the model).
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class BubbleTool : public Tool<TSeq> {
private:

    std::string _event_name;                  ///< Name of the owning intervention.
    Bubbles<TSeq> * _policy = nullptr;        ///< Resolved lazily, per model.

public:

    BubbleTool(std::string name, std::string event_name);

    /**
     * @brief Reduction applied to an exposure of this tool's agent.
     *
     * `0.0` when the transmitter shares the agent's bubble (contacts inside the
     * bubble are what the policy preserves), and `1 - f` otherwise, where `f`
     * is the intervention's transmission factor. Outside the policy window, or
     * before a partition exists, the reduction is `0.0`.
     */
    epiworld_double get_susceptibility_reduction(
        VirusPtr<TSeq> & v,
        Model<TSeq> * model
    ) override;

    std::unique_ptr<Tool<TSeq>> clone_ptr() const override;

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
 * `BubbleTool` ("Social bubble") to every agent. When a susceptible agent `p`
 * is exposed to an infectious neighbor, the tool identifies the transmitter
 * through the virus (`v->get_agent()`) and compares the two agents' bubble
 * labels:
 *
 * - same bubble  -> reduction `0.0`: contacts *inside* the bubble are what the
 *                   policy preserves, so they are left untouched;
 * - different bubbles -> reduction `1 - f`, i.e. transmission along contacts
 *                   *outside* the bubble is scaled by the transmission factor
 *                   `f`;
 * - outside the policy window -> reduction `0.0`, no effect.
 *
 * The transmission factor `f` is how strictly the bubble is observed:
 * `f == 0` is a perfectly efficient bubble (out-of-bubble contact is cut
 * entirely, which is equivalent to deleting those edges), `f == 1` disables the
 * intervention (out-of-bubble contact is as good as before), and intermediate
 * values model a soft contact reduction -- people still meet outside their
 * bubble, just less often or more carefully.
 *
 * `f` is **not** stored in the intervention: `deploy()` registers it as a model
 * parameter (`param_name`, "Bubble transmission factor" by default) and the
 * tool reads it from the model on every exposure. It can therefore be inspected
 * with `model.get_param()`, changed mid-run with `model.set_param()`, read from
 * a parameter file with `model.read_params()`, or swept over in a calibration
 * without rebuilding the intervention. Values outside `[0, 1]` are clamped.
 *
 * Since reductions combine as `1 - prod(1 - r_i)`, a reduction of `1.0`
 * (`f == 0`) zeroes the transmission probability regardless of any other tools
 * the agent carries. Contact weights are uniform in these models, so
 * suppressing transmission on out-of-bubble contacts is equivalent to deleting
 * those contacts, while keeping the network intact for other purposes (contact
 * tracing, output).
 *
 * ## Where the state lives
 *
 * `Bubbles` *is* the model's global event: `deploy()` hands the model a clone
 * of it, and the bubble partition (`get_bubble_id()`) is a member of that
 * clone. Since `Model`'s copy constructor deep-copies global events and tools,
 * every copy of a model -- including the per-thread copies `run_multiple()`
 * makes -- owns its partition outright, and the tool resolves to the
 * intervention of whichever model it belongs to. Nothing is shared between
 * models, so replicates may run on as many threads as you like.
 *
 * The object you construct is a template: after `deploy()` it is no longer
 * connected to the model, and its own `get_bubble_id()` stays empty. To look at
 * the partition of a model that has run, ask the model for its copy:
 *
 * ```cpp
 * const auto & bubble_id = Bubbles<>::get_from(model)->get_bubble_id();
 * ```
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
 * socialises with.
 *
 * ## The algorithms
 *
 * Both rules start from the **household contact graph**: one node per
 * household, with an edge between two households whenever at least one member
 * of the first is connected to a member of the second in the agents' contact
 * network. All random draws use the model's RNG, so a run is reproducible from
 * its seed.
 *
 * **`BubbleFlavor::Household`** -- grow bubbles from seed households:
 *
 * 1. Visit households in random order.
 * 2. Skip a household if it already belongs to a bubble; otherwise open a new
 *    bubble containing it, and set the *frontier* to its unassigned neighbours
 *    in the household contact graph.
 * 3. While the bubble holds fewer than `group_size` households and the frontier
 *    is not empty, draw a household from the frontier at random, add it to the
 *    bubble, and extend the frontier with that household's unassigned
 *    neighbours. A household that is tied to several members of the bubble
 *    appears in the frontier more than once and is correspondingly more likely
 *    to be drawn, so stronger ties are favoured.
 * 4. Stop when no unassigned neighbour remains, even if the bubble is smaller
 *    than `group_size`.
 *
 * Each bubble is therefore a *connected* subgraph of the household contact
 * graph -- not necessarily a clique, so with `group_size > 2` two households in
 * one bubble need not be tied to each other directly. Because a household with
 * no available partner is left on its own, the number of bubbles is at least
 * `ceil(n_households / group_size)`. Cost is linear in the number of edges.
 *
 * **`BubbleFlavor::Peer`** -- agents choose peers from those still available:
 *
 * Households are kept in a disjoint-set (union-find) structure that also tracks
 * how many households each bubble holds, so a bubble that is full can be
 * recognised at once.
 *
 * 1. Visit the agents in random order.
 * 2. Skip an agent whose household is already in a full bubble -- it has left
 *    the pool and can neither choose nor be chosen.
 * 3. Otherwise, repeatedly draw one of the agent's contacts outside its own
 *    household, at random and without replacement, until the agent has made
 *    `group_size` successful choices or has no contact left that its bubble can
 *    still take in. A draw is accepted when the two households are in different
 *    bubbles **and** the merged bubble would hold at most `max_households`
 *    households; otherwise that contact is simply unavailable and the agent
 *    draws again. Accepting a draw merges the two households (household
 *    commitment), and the agent stops once its bubble is full.
 *
 * The cap is what makes the rule work. Without it the merges percolate: with a
 * few members per household each choosing someone, the household graph becomes
 * connected and every household lands in one giant bubble -- no restriction at
 * all.
 *
 * Note that `max_households` is the effective dial on bubble size, while
 * `group_size` rarely binds: because a choice by *any* member commits the whole
 * household, the members of one household together tend to fill its bubble
 * regardless of how many choices each of them is allowed. With
 * `max_households == 2`, in particular, one accepted choice fills the bubble and
 * `group_size` has no effect at all. Households whose every contact was taken
 * first remain on their own.
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
 * applied by the intervention itself (a daily global event) and take effect the
 * following step.
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
 * // Two households per bubble from day 10, halving out-of-bubble transmission.
 * Bubbles<> bubbles(
 *     household_id, BubbleFlavor::Household,
 *     2,      // group_size
 *     0.5,    // transmission_factor (initial value of the model parameter)
 *     10      // start_day
 * );
 * bubbles.deploy(model);
 *
 * // The factor lives in the model, so it can be changed without touching
 * // the intervention.
 * model.set_param("Bubble transmission factor", 0.25);
 *
 * model.run(100, 1231);
 *
 * // The partition belongs to the model, not to `bubbles`.
 * const auto & bubble_id = Bubbles<>::get_from(model)->get_bubble_id();
 * ```
 *
 * @note Both rules produce *exclusive* bubbles, which is what the modelled
 * policies prescribe. A rule that instead grants each person a personal budget
 * of contacts that need be neither mutual nor exclusive (e.g. "up to ten
 * different people a week") is not a partition of the population and cannot be
 * expressed this way.
 *
 * @note Not yet modelled: fixed out-of-bubble "sport partners", travel between
 * regions, and caps on the number of individuals (as opposed to households) a
 * bubble may contain.
 *
 * @tparam TSeq Sequence type (should match `TSeq` across the model).
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class Bubbles final : public GlobalEvent<TSeq> {
private:

    std::vector< size_t > household_id;
    BubbleFlavor flavor;
    size_t group_size;               ///< households per bubble (Household) or max peers (Peer).
    size_t max_households;           ///< cap on households per bubble (Peer only).
    epiworld_double transmission_factor; ///< initial value of the model parameter, in [0, 1].
    int start_day;
    int end_day;
    int rewire_every;
    std::string param_name;          ///< model parameter holding the transmission factor.

    // Per-model state. Copied with the intervention, so each model (including
    // each per-thread copy made by run_multiple) owns its own partition.
    std::vector< int > bubble_id;    ///< Per-agent bubble label; -1 = unassigned.
    int last_sim_id = -1;            ///< Sim id the current partition was computed for.
    int last_epoch  = -1;            ///< Rewiring epoch the current partition was computed for.

    void partition_household(Model<TSeq> * model);
    void partition_peer(Model<TSeq> * model);

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
     * @param transmission_factor Initial value of the model parameter
     *        `param_name`: the multiplier applied to transmission between
     *        agents of *different* bubbles, in `[0, 1]`. `0.0` (the default) is
     *        a perfectly efficient bubble -- contact outside it is cut
     *        entirely; `0.5` halves out-of-bubble transmission (a soft contact
     *        reduction); and `1.0` turns the intervention off. Transmission
     *        within a bubble is never altered. The value is stored in the model
     *        (see `deploy()`), not here, so it can be changed at any time with
     *        `model.set_param(param_name, ...)`.
     * @param start_day First day on which the policy applies.
     * @param end_day Day on which the policy is lifted (exclusive). Use a
     *        negative value for a policy that never ends. Must be greater than
     *        `start_day`.
     * @param rewire_every Re-randomise the bubbles every this-many days, for
     *        policies whose contacts change over time. `0` keeps the bubbles
     *        fixed for the whole intervention.
     * @param name Name given to the tool and to the intervention's global
     *        event; the tool finds its policy by this name, and `get_from()`
     *        looks it up on a model with it.
     * @param max_households **`Peer` flavor only**: the largest number of
     *        households a bubble may contain. A nomination that would exceed it
     *        is declined, which is what stops one household's choices from
     *        chaining into the next (see `BubbleFlavor::Peer`). The default of
     *        `2` represents two households joined by a close contact. Ignored by
     *        the `Household` flavor, where `group_size` already is the cap.
     * @param param_name Name of the model parameter that holds the transmission
     *        factor. Give two interventions deployed on the same model
     *        different names if they are to be dialled independently.
     *
     * @throws std::range_error if `transmission_factor` is outside `[0, 1]`, if
     *         `group_size` is zero for the `Household` flavor, if
     *         `max_households` is less than 2 for the `Peer` flavor, or if
     *         `end_day` is non-negative and not greater than `start_day`.
     */
    Bubbles(
        std::vector< size_t > household_id,
        BubbleFlavor flavor,
        size_t group_size,
        epiworld_double transmission_factor = 0.0,
        int start_day = 0,
        int end_day = -1,
        int rewire_every = 0,
        std::string name = "Social bubble",
        size_t max_households = 2u,
        std::string param_name = "Bubble transmission factor"
    );

    /**
     * @brief Install the intervention on a model.
     *
     * Registers the model parameter `param_name` (set to the
     * `transmission_factor` passed to the constructor, overwriting any value it
     * already had), hands the model a copy of this intervention as a global
     * event, and adds the bubble `Tool`, which is given to every agent at reset
     * time -- the same moment the partition is computed.
     *
     * To drive the factor from a parameter file instead, call
     * `model.read_params()` (or `set_param()`) *after* `deploy()`.
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
     * @brief The model's copy of a deployed intervention.
     *
     * @param model Model the intervention was deployed on.
     * @param name Name it was deployed under.
     * @return Pointer to the model's own intervention, or `nullptr` if the
     *         model has no global event by that name, or it is not a `Bubbles`.
     */
    static Bubbles<TSeq> * get_from(
        Model<TSeq> & model,
        const std::string & name = "Social bubble"
    );

    /**
     * @brief Recompute the partition using the model's RNG.
     *
     * Called at reset time (through the tool's distribution function) and at
     * each rewiring epoch. Rarely needed directly.
     */
    void compute_partition(Model<TSeq> * model);

    /**
     * @brief Current bubble label of every agent, indexed by agent id.
     *
     * Two agents in the same bubble transmit as they would without the policy;
     * transmission between labels is scaled by the transmission factor.
     * Populated on the first reset, and only on the model's own copy of the
     * intervention -- see `get_from()`.
     */
    const std::vector< int > & get_bubble_id() const;

    /// @brief The rule used to form bubbles.
    BubbleFlavor get_flavor() const;

    /**
     * @brief Name of the model parameter holding the transmission factor.
     *
     * Use it with `model.get_param()` / `model.set_param()` to read or change
     * how leaky the bubbles are, including in the middle of a run.
     */
    const std::string & get_param_name() const;

    /// @brief True when the policy applies on the given day.
    bool is_active(int today) const;

    /**
     * @brief Rewiring epoch the current partition was computed for.
     *
     * `0` for the partition drawn at reset, incrementing every `rewire_every`
     * days while the policy is active. `-1` before the first partition.
     */
    int get_last_epoch() const;

    /**
     * @brief Reduction applied to an exposure between two agents.
     *
     * The tool's side of the intervention; see `BubbleTool`.
     */
    epiworld_double susceptibility_reduction(
        const Agent<TSeq> * p,
        const Agent<TSeq> * transmitter,
        Model<TSeq> * model
    ) const;

    /// @brief Re-randomises the partition at rewiring epochs.
    void operator()(Model<TSeq> * model, int day) override;

    std::unique_ptr< GlobalEvent<TSeq> > clone_ptr() const override;

};

#endif
