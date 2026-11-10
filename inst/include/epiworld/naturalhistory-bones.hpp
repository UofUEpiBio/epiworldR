#ifndef EPIWORLD_NATURALHISTORY_BONES_HPP
#define EPIWORLD_NATURALHISTORY_BONES_HPP

/**
 * @file naturalhistory-bones.hpp
 * @brief The natural history of a disease: its stages and the declared
 * transitions between them (see `natural-history.md`).
 *
 * @details A natural history is a small graph. Its nodes are **stages**
 * (susceptible, latent, infectious, recovered, ...) and its edges are
 * **transitions**, each of one kind:
 *
 * - `infection()`: by contact, drawn from the model's layers (or a custom
 *   sampler). Taking it gives the agent the virus.
 * - `rate()`: a daily probability. The rates out of a stage compete: one
 *   draw (`roulette()`) picks one of them, or none.
 * - `after()`: a duration, drawn when the agent enters the stage. It is a
 *   deadline: the agent leaves on the update of the day it falls, unless a
 *   rate fired first.
 * - `external()`: fired by code, through `Model::fire()`.
 *
 * A stage is a **carrier** (its agents hold the virus) and may be
 * **infectious** (they transmit; this implies carrier). Entering a
 * non-carrier stage releases the virus. There are no fixed roles: a stage is
 * susceptible because an infection transition leaves it, and a death is a
 * stage with contact weight 0 and no exits.
 *
 * Every number a transition uses is a `Value`: a number, a model parameter
 * (by name) or a function. A number given directly becomes a model parameter
 * when the natural history is set on a model.
 */

class Nodes;

class NaturalHistory;

class TransitionRef;

/// Flags of a stage (see `NaturalHistory::add_stage()`).
namespace stage {
    constexpr unsigned int carrier    = 1u;      ///< Agents in it hold the virus
    constexpr unsigned int infectious = 1u | 2u; ///< They transmit (implies carrier)
}

/// Which tool effect a rate takes (see `TransitionRef::effect()`).
enum class Effect : uint8_t {
    none,     ///< No tool effect
    recovery, ///< `p` becomes `1 - (1 - p) * (1 - recovery enhancer)`
    death     ///< `p` becomes `p * (1 - death reduction)`
};

/// One or more stages (an `int` or a list of them).
class Nodes {
public:
    std::vector< int > ids;
    Nodes(int id) : ids({id}) {}
    Nodes(std::initializer_list< int > l) : ids(l) {}
    Nodes(std::vector< int > v) : ids(std::move(v)) {}
};

using ValueFun = std::function< epiworld_double(Agent &, Model &) >;

using DurationFun = std::function< int(Agent &, Model &) >;

/// A custom infection sampler: the virus the agent catches, or `nullptr`.
using InfectionSampler = std::function< Virus * (Agent &, Model &) >;

/**
 * @brief A number a transition uses: a constant, a model parameter, or a
 * function.
 * @details
 * - `Value(0.1)`: a number. When the natural history is set on a model, it
 *   becomes a model parameter (named after the transition).
 * - `Value("Recovery rate")`: a model parameter, by name.
 * - `Value::inverse("Days")`, `Value::complement("p")`: `1 / p`, `1 - p`.
 * - `Value(fun)`: `fun(agent, model)`.
 *
 * A parameter is read through the agent's virus, so a variant can point it to
 * another parameter (`Virus::use_param()`).
 */
class Value {
public:
    enum class Kind : uint8_t { none, number, param, inverse, complement, fun };

private:
    Kind kind = Kind::none;
    epiworld_double number = 0.0;
    ParamRef param = ParamRef("");
    ValueFun fun = nullptr;

    friend class NaturalHistory;

public:
    Value() = default;
    /// A number (any arithmetic type, whatever `epiworld_double` is)
    template<
        typename T,
        typename std::enable_if< std::is_arithmetic< T >::value, int >::type = 0
    >
    Value(T x) : kind(Kind::number), number(static_cast< epiworld_double >(x)) {}
    Value(const char * name) : kind(Kind::param), param(std::string(name)) {}
    Value(std::string name) : kind(Kind::param), param(std::move(name)) {}
    Value(ValueFun f) : kind(Kind::fun), fun(std::move(f)) {}

    static Value inverse(std::string name);    ///< `1 / parameter`
    static Value complement(std::string name); ///< `1 - parameter`

    Kind get_kind() const { return kind; }
    bool is_set() const { return kind != Kind::none; }

    /// The parameter's name ("" for a function).
    const std::string & get_param_name() const { return param.name(); }

    /// The value for this agent (its virus may redirect the parameter).
    epiworld_double eval(Agent & agent, Model & model) const;

    /// The same, with the parameter redirected by `virus` (e.g., the
    /// infector's, for a probability of infection).
    epiworld_double eval(
        Agent & agent,
        const Virus * virus,
        Model & model
    ) const;
};

/**
 * @brief How many days a stage lasts, drawn when the agent enters it.
 * @details Always a whole number of days, at least one.
 * - `fixed(p)`: `p`, rounded.
 * - `geometric(mean)`: a geometric number of days with that mean (the same as
 *   leaving each day with probability `1 / mean`).
 * - `lognormal(mean, sd)`, `gamma(mean, sd)`: the mean and standard deviation
 *   of the duration, in days. A standard deviation of zero is fixed.
 * - `fun`: `fun(agent, model)`.
 *
 * `since(s)` counts the duration from the day the agent entered stage `s`
 * (in the same infection), instead of the day it entered this stage.
 */
class Duration {
public:
    enum class Kind : uint8_t { none, fixed, geometric, lognormal, gamma, fun };

private:
    Kind kind = Kind::none;
    Value mean;
    Value sd;
    DurationFun fun = nullptr;
    int anchor = -1;

    friend class NaturalHistory;

public:
    Duration() = default;
    Duration(DurationFun f) : kind(Kind::fun), fun(std::move(f)) {}

    static Duration fixed(Value days);
    static Duration geometric(Value mean);
    static Duration lognormal(Value mean, Value sd);
    static Duration gamma(Value mean, Value sd);

    /// Counts from the day the agent entered stage `s` in this infection.
    Duration since(int s) const;

    Kind get_kind() const { return kind; }
    int get_anchor() const { return anchor; }

    /// Draws a number of days (at least one).
    int draw(Agent & agent, Model & model) const;
};

/// The kind of a transition (see `naturalhistory-bones.hpp`).
enum class TransitionKind : uint8_t { infection, rate, after, external };

/**
 * @brief A declared transition: from a set of stages to a stage (or a branch
 * of stages), of one kind.
 */
struct Transition {
    std::string name;
    TransitionKind kind;
    std::vector< int > from;
    int to = -1;

    /// A branch: the target is drawn among these, with these probabilities.
    std::vector< int > branch_to;
    std::vector< Value > branch_prob;

    Value value;          ///< infection (per contact) and rate
    Duration duration;    ///< after
    Effect effect = Effect::none;
    InfectionSampler sampler = nullptr; ///< infection only (optional)

    /// A status transition: only for agents in one of these stages (empty: any)
    std::vector< int > when;
};

/**
 * @brief Refines a transition just declared (chainable); converts to its id.
 */
class TransitionRef {
private:
    NaturalHistory * nh;
    int id;

public:
    TransitionRef(NaturalHistory * nh_, int id_) : nh(nh_), id(id_) {}

    /// The transition's id for `Model::fire()` and reactions: a stage
    /// transition's own index, or, for a status transition, `-2 - index`.
    operator int() const;
    int get_id() const { return static_cast< int >(*this); }

    TransitionRef & name(std::string n);     ///< Default: "From -> To"
    TransitionRef & effect(Effect e);        ///< The tool effect a rate takes

    /// The target is drawn among `targets`, with their probabilities.
    TransitionRef & to(std::vector< std::pair< int, Value > > targets);

    /// A custom infection sampler (infection transitions only; C++ only).
    TransitionRef & sampler(InfectionSampler f);

    /// A status transition: only for agents in these stages (read at the
    /// start of the day).
    TransitionRef & when(Nodes stages);
};

/**
 * @brief What the engine needs about a stage, derived from the declarations.
 */
struct StageInfo {
    int infection = -1;           ///< The infection transition out of it
    std::vector< int > rates;     ///< Its rate transitions, in order
    int after = -1;               ///< Its duration transition
    bool needs_visit = false;     ///< Has a rate or a duration
    bool pushable = false;        ///< Its only drawn exit is a default infection
};

/**
 * @brief The natural history of a disease: stages and declared transitions.
 * @details It is built on the virus (`Virus::add_stage()` and the others),
 * and the model takes its stages from `Model::add_virus()`. The first stage
 * is where agents start. Example (SIR):
 *
 * @code{.cpp}
 * Virus flu("Flu");
 * int S = flu.add_stage("Susceptible");
 * int I = flu.add_stage("Infected", stage::infectious);
 * int R = flu.add_stage("Recovered");
 * flu.infection(S, I, "Transmission rate");
 * flu.rate(I, R, "Recovery rate").effect(Effect::recovery);
 * model.add_virus(flu, assign::share(0.01));
 * @endcode
 */
class NaturalHistory {
protected:
    std::string name;
    std::vector< std::string > stage_labels;
    std::vector< unsigned int > stage_flags;
    std::vector< epiworld_double > stage_contact;
    std::vector< Transition > transitions;

    // Derived by finalize()
    std::vector< StageInfo > stages;
    bool finalized = false;

    friend class TransitionRef;
    friend class Model;
    bool status_chart = false; ///< See `Statuses`

    int add_transition(Transition t);
    void check_stage(int s, const char * what) const;

    /// Checks the declarations and derives `stages`; turns numbers into model
    /// parameters, named after `owner` (the virus) and the transition.
    void finalize(Model & model, const std::string & owner);

public:
    explicit NaturalHistory(std::string name_);

    /// Adds a stage; returns its id. Flags: `stage::carrier`, `stage::infectious`.
    int add_stage(std::string label, unsigned int flags = 0u);

    /// Contact weight of a stage (default 1). 0 takes its agents out of every
    /// layer: they are not met and do not transmit.
    void set_contact(int s, epiworld_double weight);

    TransitionRef infection(Nodes from, int to, Value per_contact);
    TransitionRef rate(Nodes from, int to, Value daily_probability);
    TransitionRef after(int from, int to, Duration d);
    TransitionRef external(Nodes from, int to, std::string label);

    const std::string & get_name() const { return name; }
    size_t get_n_stages() const { return stage_labels.size(); }
    size_t get_n_transitions() const { return transitions.size(); }
    const std::vector< std::string > & get_stages() const { return stage_labels; }
    const Transition & get_transition(int id) const { return transitions.at(id); }

    bool is_carrier(int s) const { return (stage_flags.at(s) & stage::carrier) != 0u; }
    bool is_infectious(int s) const { return (stage_flags.at(s) & 2u) != 0u; }
    epiworld_double get_contact(int s) const { return stage_contact.at(s); }

    /// The stage's id, by label (throws if absent).
    int stage_id(std::string_view label) const;

    /// The transition's id, by name (throws if absent).
    int transition_id(std::string_view label) const;

    /// Whether `s` is one of the transition's sources.
    bool is_source(int transition, int s) const;

    /**
     * @brief Whether `other` has the same structure: the same stages (labels,
     * flags, contact weights) and the same transitions (names, kinds,
     * sources and targets). Their values may differ.
     */
    bool same_structure(const NaturalHistory & other) const;

    /// What the engine derived about stage `s` (once added to a model).
    const StageInfo & get_stage_info(int s) const { return stages.at(s); }

    /**
     * @brief The stages and transitions as the body of a mermaid flowchart.
     * @details One node per stage (id `prefix` + index; infectious, carrier
     * and contact-0 stages get the classes `infectious`, `carrier` and
     * `away`) and one unlabeled arrow per pair of nodes that some transition
     * joins. The kinds and values of the transitions are left out: they are
     * the model's parameters, and differ between models built on the same
     * stages. See `Model::draw_natural_history()`.
     * @param prefix Prefix of the node ids.
     */
    std::string draw_mermaid(const std::string & prefix) const;
};

/**
 * @brief The statuses of a model: where agents are (free, quarantined,
 * isolated, hospitalized), and the transitions between them.
 * @details The same vocabulary as a natural history (`rate()`, `after()`,
 * `external()`), without infections or carriers. A status transition can be
 * restricted to agents in some stages (`TransitionRef::when()`); a duration
 * can count from a stage's entry (`Duration::since()`). A status with contact
 * weight 0 takes its agents out of every layer. The first status is where
 * agents start. Set on a model with `Model::set_statuses()`; a model without
 * statuses has one, "Free".
 */
class Statuses : public NaturalHistory {
public:
    Statuses() : NaturalHistory("Statuses") { this->status_chart = true; }

    /// Adds a status; returns its id.
    int add_status(std::string label) { return this->add_stage(std::move(label)); }

    TransitionRef infection(Nodes, int, Value) = delete;
};

/**
 * @brief What a reaction responds to (see `Model::react()`).
 * @details `Trigger::entered(s)`: an agent enters stage `s`;
 * `Trigger::fired(t)`: transition `t` is applied. `.was({statuses})`: only if
 * the agent's status at the start of the day was one of these.
 */
struct Trigger {
    enum class Kind : uint8_t { entered, fired } kind = Kind::entered;
    int target = -1;
    std::vector< int > was_status; ///< Empty: any status

    static Trigger entered(int stage) { return {Kind::entered, stage, {}}; }
    static Trigger fired(int transition) { return {Kind::fired, transition, {}}; }
    Trigger was(std::vector< int > statuses) const
    {
        Trigger t(*this);
        t.was_status = std::move(statuses);
        return t;
    }
};

#endif
