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

/**
 * @brief How a bubble is *realized* in the contact network.
 * @ingroup globalevents
 *
 * @details Orthogonal to `BubbleFlavor`, which decides how households are
 * grouped. This decides what being in a bubble together actually does.
 */
enum class BubbleTies {
    /**
     * The contact network is never modified. Being in a bubble together only
     * means transmission between the two is left alone, while contact outside
     * the bubble is damped by the transmission factor.
     *
     * This makes the policy purely subtractive, which under-states it: two
     * households merged because one member of each happened to be tied keep
     * exactly that one contact between them, and everybody else in the two
     * households never meets. The same gap exists inside a household, where two
     * members with no edge between them never meet either.
     */
    Existing,
    /**
     * The bubble is completed to a clique: every pair of agents sharing a bubble
     * is tied for as long as the bubble lasts, and the ties are withdrawn when
     * it is redrawn or the policy lifts.
     *
     * This is what the modelled policy actually says -- the household you have
     * bubbled with is a household you now see. Two consequences worth knowing:
     * households are completed too (members who were not tied in the underlying
     * network now are), and the force of infection inside a bubble grows with
     * bubble size, since every tie is an independent daily exposure.
     *
     * Contacts *outside* the bubble are untouched and still damped by the
     * transmission factor, which is what makes the bubble imperfect.
     *
     * The ties only exist while the policy is in force. If household members
     * should meet before and after it too, those ties belong in the contact
     * network itself. For the same reason `Complete` is not neutral at a
     * transmission factor of `1`: it adds contact inside the bubble and damps
     * nothing outside it, so transmission rises above running without the
     * policy.
     */
    Complete
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
 * Users do not instantiate this directly; the intervention hands it to every
 * agent when it sets itself up, at the start of each run.
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
 * The intervention is a global event, so installing it is a one-liner:
 *
 * ```cpp
 * model.add_globalevent(bubbles);
 * ```
 *
 * Everything else happens by itself. At the start of every run -- from
 * `Model::reset()`, so the policy is in force on day 1 -- the intervention sets
 * itself up on the model it is running in: it registers the
 * transmission-factor parameter, draws the bubble partition with that run's
 * RNG, and hands a `BubbleTool` ("Social bubble") to every agent. Each
 * replicate of `run_multiple()` -- each of which runs on its own copy of the
 * model, possibly on its own thread -- repeats this with its own seed.
 *
 * By default the contact network is **not modified** (`BubbleTies::Existing`);
 * with `BubbleTies::Complete` the bubble is additionally realized as ties, which
 * is described under *Two ways to be in a bubble* below. Either way, when a
 * susceptible agent `p` is exposed to an infectious neighbor, the tool
 * identifies the transmitter through the virus (`v->get_agent()`) and compares
 * the two agents' bubble labels:
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
 * `f` is **not** stored in the intervention: setup registers it as a model
 * parameter (`param_name`, "Bubble transmission factor" by default) and the
 * tool reads it from the model on every exposure. It can therefore be inspected
 * with `model.get_param()`, changed mid-run with `model.set_param()`, or swept
 * over in a calibration without rebuilding the intervention. The value given to
 * the constructor is only a default: if the model already carries that
 * parameter -- because it was set with `model.add_param()` or read from a
 * parameter file with `model.read_params()` before the run -- the model's value
 * is kept. Values outside `[0, 1]` are clamped.
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
 * `Bubbles` *is* the model's global event: `add_globalevent()` hands the model a
 * clone of it, and the bubble partition (`get_bubble_id()`) is a member of that
 * clone. Since `Model`'s copy constructor deep-copies global events and tools,
 * every copy of a model -- including the per-thread copies `run_multiple()`
 * makes -- owns its partition outright, and the tool resolves to the
 * intervention of whichever model it belongs to. Nothing is shared between
 * models, so replicates may run on as many threads as you like.
 *
 * The object you construct is a template: once added to a model it is no longer
 * connected to it, and its own `get_bubble_id()` stays empty. To look at the
 * partition of a model that has run, ask the model for its copy:
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
 * The last point is essential, and it holds for both realizations. Under
 * `BubbleTies::Existing` the intervention can only suppress transmission along
 * edges that exist, so putting two households that share no contact into one
 * bubble would change nothing at all and `group_size` would be inert. Under
 * `BubbleTies::Complete` the ties would be created, but bubbling two households
 * that have never met is not the policy being modelled: a household chooses a
 * bubble partner it already socialises with.
 *
 * ## Two ways to be in a bubble
 *
 * `BubbleFlavor` decides *how households are grouped*. `BubbleTies` decides
 * *what being grouped does*, and the two are independent.
 *
 * Under **`BubbleTies::Existing`** (the default, and the original behaviour) a
 * bubble is only a transmission rule. This under-states the policy in a way
 * worth being explicit about: if households `a` and `b` were merged because one
 * member of each happened to be tied, that single contact is the only one
 * preserved between them -- the rest of `a` never meets the rest of `b`, though
 * the policy says the two households have bubbled. The same gap exists inside a
 * household, whose members need not all be tied to each other.
 *
 * Under **`BubbleTies::Complete`** the bubble is completed to a clique: every
 * pair of agents sharing a bubble is tied for as long as the bubble lasts.
 * Households become complete, and merged households meet in full. Ties the
 * network already had are left alone; the intervention records only the ones it
 * created, and withdraws exactly those when the bubble is redrawn, when the
 * policy lifts, and on the last day of the run -- so a run never leaves the
 * network changed. `get_created_ties()` lists them while they are up, and
 * `restore_network()` withdraws them early if a hand-driven day loop stops
 * before the end.
 *
 * Two consequences to keep in mind. Completing a bubble **raises the force of
 * infection inside it**: every tie is an independent daily exposure, so a bubble
 * of `k` agents gives each member `k - 1` chances a day rather than whatever
 * degree they had. And the ties are **real while they are up** -- they show up
 * in `write_edgelist()`, in `get_n_neighbors()`, and in contact tracing.
 *
 * Contacts *outside* the bubble are untouched by either setting: they stay in
 * the network and are damped by `f`, which is what makes the bubble imperfect.
 *
 * The ties only exist while the policy is in force: households are completed
 * on `start_day` and go back to however the network ties them when it lifts.
 * If household members should meet before and after the policy too, put those
 * ties in the contact network itself. For the same reason `Complete` is not
 * neutral at `f == 1`: it adds contact inside the bubble and damps nothing
 * outside it, so transmission rises above running without the policy.
 *
 * Grouping never sees these ties. Both rules read the contact network to decide
 * which households may bubble together, and they skip every tie a bubble policy
 * on the model is holding -- this one's from the previous epoch, or another
 * policy's -- so `BubbleFlavor` draws the same bubbles whatever `BubbleTies`
 * says, and one policy's grouping does not depend on another's ties.
 *
 * Two `Complete` policies may run on the same model -- household bubbles and
 * school bubbles, say. A tie both want is on the books of whichever created it,
 * and when that policy's window closes the tie is handed to the other rather
 * than removed, so a bubble that is still open does not lose a contact because
 * a different policy ended.
 *
 * Anything else that edits the network is left to it. If another event takes
 * away a tie the bubble wants -- isolating an agent, say -- the tie stays gone
 * for as long as that bubble stands: the policy does not put it back. A tie the
 * policy created stays on its books all the same, so if it is restored while
 * the bubble is up (the isolation ends) it still comes down with the bubble.
 * The policy only ever takes back what it put in, and never puts back what
 * something else took out. A bubble redrawn at a later epoch is a new bubble,
 * though, and is completed afresh; an event that means to keep agents apart
 * across a redraw has to say so again.
 *
 * `BubbleTies::Complete` needs an undirected model and cannot be combined with
 * `Model::set_rewire_fun()` -- a rewiring function moves ties between agents,
 * so a tie the intervention created could be moved out from under it and never
 * withdrawn. Both are checked when the intervention sets itself up. The
 * rewiring check is on the function, not on the proportion: `Model::rewire()`
 * calls the function on every step.
 *
 * Finally, the clique has to fit the virus sampler, which weighs at most half of
 * `Model::array_double_tmp` in neighbors at once when pulling (`roulette()` uses
 * two slots per candidate). The policy answers for the ties it adds, not for the
 * network it was given: a bubble is refused if completing it would take a
 * member past that ceiling, or add to one already past it, while a member the
 * network itself put past it -- a hub -- is fine as long as its bubble adds
 * nothing to it. What binds is each member's resulting *degree*, not the size of
 * the bubble: the clique is added on top of the ties an agent already has
 * outside it.
 *
 * ## The algorithms
 *
 * Both rules start from the **household contact graph**: one node per
 * household, with an edge between two households whenever at least one member
 * of the first is connected to a member of the second in the agents' contact
 * network (ties a `Complete` bubble is holding do not count; see above). All
 * random draws use the model's RNG, so a run is reproducible from its seed.
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
 * The initial partition is drawn when the intervention sets itself up, at reset
 * time, so it is redrawn for each replicate of `run_multiple()` using that
 * replicate's seed and is already in force on day 1. Re-randomisations are
 * applied by the intervention itself (a daily global event) and, since global
 * events run after each day's transitions, take effect the following step.
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
 *     0.5,    // transmission_factor (default of the model parameter)
 *     10      // start_day
 * );
 * model.add_globalevent(bubbles);
 *
 * // The factor lives in the model, so it can be changed without touching the
 * // intervention. Setting it before the run overrides the default above.
 * model.add_param(0.25, "Bubble transmission factor", true);
 *
 * model.run(100, 1231);
 *
 * // The partition belongs to the model, not to `bubbles`.
 * const auto & bubble_id = Bubbles<>::get_from(model)->get_bubble_id();
 * ```
 *
 * To have those two households actually meet -- every member of one tied to
 * every member of the other, and each household complete in itself -- ask for
 * the bubble to be realized as ties:
 *
 * ```cpp
 * bubbles.set_ties(BubbleTies::Complete);
 * model.add_globalevent(bubbles);
 * ```
 *
 * A household-only lockdown in which the household is complete, and everyone
 * else is reachable only through an imperfectly observed bubble, is
 * `group_size == 1` with a non-zero factor:
 *
 * ```cpp
 * Bubbles<> household_only(
 *     household_id, BubbleFlavor::Household,
 *     1,      // group_size: nobody bubbles with another household
 *     0.15,   // but the bubble leaks: 15% of out-of-bubble transmission remains
 *     0, -1, 0, "Social bubble", 2u, "Bubble transmission factor",
 *     BubbleTies::Complete
 * );
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
    BubbleTies ties;                 ///< whether the bubble is completed to a clique.

    // Per-model state. Copied with the intervention, so each model (including
    // each per-thread copy made by run_multiple) owns its own partition.
    std::vector< int > bubble_id;    ///< Per-agent bubble label; -1 = unassigned.

    /**
     * @brief The ties this intervention created, and only those.
     *
     * @details Completing a bubble mostly re-uses ties the network already had;
     * `Model::add_edge()` reports which ones were genuinely new, and only those
     * are recorded here. Withdrawing a bubble then removes exactly what it
     * added, leaving anything the model itself did to the network alone -- which
     * is why this is a list of ties rather than, say, a remembered degree.
     */
    std::vector< std::pair< size_t, size_t > > created_ties;
    int ties_epoch = -1;             ///< Epoch the materialized ties belong to.
    int model_id   = -1;             ///< Sim id this intervention was set up for.
    int last_epoch = -1;             ///< Rewiring epoch the current partition was computed for.

    /**
     * @brief The ties bubble policies on a model are holding -- this one's and
     * any other's -- looked up one agent at a time.
     *
     * @details `focus(a)` marks the held partners of agent `a`, after which
     * `contains(b)` says whether the tie `a`--`b` is held. Both cost no more
     * than the ties involved, so a grouping rule can walk every contact of
     * every agent and skip the held ones for the price of the walk, however
     * many ties are held. With nothing held, both are no-ops.
     */
    class HeldTies {
    public:
        HeldTies(Model<TSeq> * model, const Bubbles<TSeq> & self);
        void focus(size_t a);
        bool contains(size_t b) const;
    private:
        std::vector< size_t > start;    ///< Offsets into `partner`, per agent.
        std::vector< size_t > partner;  ///< Held partners, agent by agent.
        std::vector< char > marked;     ///< Partners of the focused agent.
        size_t focused = 0u;
    };

    /**
     * @name The two grouping rules.
     *
     * @param held Ties bubble policies are holding; they are not contacts, so
     *        the rules skip them (see `compute_partition()`).
     */
    ///@{
    void partition_household(Model<TSeq> * model, HeldTies & held);
    void partition_peer(Model<TSeq> * model, HeldTies & held);
    ///@}

    /**
     * @brief Install the intervention on the model it is running in.
     *
     * @details Called from `reset()`, i.e. once per run and once per replicate
     * of `run_multiple()`. It registers the transmission-factor parameter
     * (keeping any value the model already has), draws the epoch-0 partition
     * with the run's RNG, and gives the bubble tool to every agent.
     *
     * @param model Model the intervention belongs to.
     * @throws std::length_error if `household_id` does not have exactly one
     *         entry per agent.
     */
    void _setup(Model<TSeq> * model);

    /// @brief Ties every pair sharing a bubble, recording what it created.
    void build_ties(Model<TSeq> * model);

    /**
     * @brief Makes the network match the policy as of the next simulation step.
     *
     * @details Builds the clique when the policy is about to apply and it is not
     * up already, and drops it when the policy is about to lapse. At a
     * rewiring epoch the daily event withdraws the old clique before drawing
     * the new partition, so this builds the new one. A standing clique is not
     * re-checked: a tie something else took away stays gone (though it stays
     * on the books).
     */
    void sync_ties(Model<TSeq> * model);

    /// @brief Whether the clique should be up for the next simulation step.
    bool wants_ties(Model<TSeq> * model) const;

    /**
     * @brief Withdraw the ties this policy is holding.
     *
     * @details A tie that is no longer in the network -- something else took
     * it away -- is not an error: it is handed over like any other, or skipped
     * if nobody wants it.
     *
     * @param hand_over When true, a tie that another `Bubbles` policy still
     *        wants is transferred to it rather than removed, so a bubble that is
     *        still open does not lose a contact because a different policy's
     *        window closed.
     */
    void withdraw_ties(Model<TSeq> * model, bool hand_over);

    /// @brief Another active `Complete` policy that wants the tie `i`--`j`.
    Bubbles<TSeq> * heir_of(Model<TSeq> * model, size_t i, size_t j);

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
     * @param ties Whether the bubble is only a transmission rule
     *        (`BubbleTies::Existing`, the default) or is completed to a clique
     *        of temporary ties (`BubbleTies::Complete`). See `BubbleTies`.
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
        std::string param_name = "Bubble transmission factor",
        BubbleTies ties = BubbleTies::Existing
    );

    /**
     * @brief The model's copy of the intervention.
     *
     * @param model Model the intervention was added to.
     * @param name Name it was added under.
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
     * Called when the intervention sets itself up (at reset) and at each
     * rewiring epoch. Rarely needed directly.
     */
    void compute_partition(Model<TSeq> * model);

    /**
     * @brief Current bubble label of every agent, indexed by agent id.
     *
     * Two agents in the same bubble transmit as they would without the policy;
     * transmission between labels is scaled by the transmission factor.
     * Populated at the start of each run, and only on the model's own copy of
     * the intervention -- see `get_from()`.
     */
    const std::vector< int > & get_bubble_id() const;

    /// @brief The rule used to form bubbles.
    BubbleFlavor get_flavor() const;

    /// @brief Whether the bubble is completed to a clique of temporary ties.
    BubbleTies get_ties() const;

    /**
     * @brief Choose how the bubble is realized. See `BubbleTies`.
     *
     * Set it before the run; the intervention reads it when it sets itself up.
     * Remember that the model owns its own copy once the intervention has been
     * added, so set this on the object you are about to add, or on the model's
     * copy via `get_from()`.
     */
    void set_ties(BubbleTies ties);

    /**
     * @brief Ties this intervention is currently keeping up, if any.
     *
     * Empty under `BubbleTies::Existing`, and empty under `Complete` outside
     * the policy window. Each pair is a tie a bubble policy created -- this
     * one, or another that handed it over when its own window closed. Ties the
     * network already had are not listed, because they are not the
     * intervention's to withdraw.
     */
    const std::vector< std::pair< size_t, size_t > > & get_created_ties() const;

    /**
     * @brief Withdraw the ties the intervention is keeping up.
     *
     * @details A run withdraws them by itself on its last day, so this is only
     * needed when the day loop is driven by hand and stops while the policy is
     * still in force. The next `Model::reset()` withdraws them too; call this
     * to get the original network back sooner. Harmless if there is nothing to
     * withdraw.
     */
    void restore_network(Model<TSeq> * model);

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
     * `0` for the partition drawn at setup, incrementing every `rewire_every`
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

    /**
     * @brief Installs the intervention on the model, once per run.
     *
     * Called by `Model::reset()`; see `_setup()`. This is what makes adding the
     * intervention to a model the only step there is.
     */
    void reset(Model<TSeq> * model) override;

    /// @brief Re-randomises the partition at rewiring epochs, and keeps the
    /// ties of a `Complete` bubble in step with the policy.
    void operator()(Model<TSeq> * model, int day) override;

    std::unique_ptr< GlobalEvent<TSeq> > clone_ptr() const override;

};

#endif
