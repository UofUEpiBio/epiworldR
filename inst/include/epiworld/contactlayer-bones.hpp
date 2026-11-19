#ifndef EPIWORLD_CONTACTLAYER_BONES_HPP
#define EPIWORLD_CONTACTLAYER_BONES_HPP

/**
 * @file contactlayer-bones.hpp
 * @brief Base class of the contact layers of a model.
 *
 * @details A *contact layer* is one source of contacts of a model: it answers
 * "whom did agent `i` contact this step?". A model owns any number of them (see
 * `Model::add_layer()`), and none are implied: a model without layers has no
 * transmission through contacts.
 *
 * The model drives the life cycle of its layers:
 *
 * - `reset()` at the start of `Model::reset()`, before anything is distributed
 *   (viruses, tools, the global events' set-up), so that those see
 *   the run's contacts. A layer restores its baseline if it was edited during
 *   the previous run (see `is_dirty()`), or draws a new one, and sizes its
 *   per-run buffers.
 * - `update()` at the end of every step, including day 0 (at the start of
 *   `Model::next()`): after the agents' update functions, the global events
 *   and rewiring, so that whatever changed during a step (states, edits of the
 *   layer) is in force for the next step's contacts.
 * - `clone()` whenever the model is copied (`run_multiple()` copies models),
 *   so replicates never share layer state.
 *
 * ### Ownership
 *
 * Layers are created as shared pointers (`layer::network()`,
 * `layer::mixing()`; `std::make_shared` for a layer type of your own), and
 * `Model::add_layer()` shares the pointer: the model and whoever else holds it
 * see the same layer, so an edit made through either is an edit of the model.
 * A layer belongs to one model at most. A copy of the model clones its layers
 * (the pointers held outside keep pointing at the original's), and a layer
 * outlives its model only detached from it. A copy of a layer is detached too.
 *
 * ### Weight and viruses
 *
 * `set_weight()` makes the layer's contacts count more or less than the others':
 * the odds of infection through the layer are multiplied by the weight (a
 * model parameter, read at every step; 1 if none is set; 0 turns the layer
 * off). `add_target()` restricts the layer to some diseases, by name, as tools
 * are restricted; by default every disease spreads through every layer.
 *
 * ### Edits and the dirty flag
 *
 * Whether an edit is permanent depends only on when it is made. An edit made
 * during a run (from an update function or a global event) applies to that run
 * only: it sets the *dirty* flag, and `reset()` restores the baseline of a
 * dirty layer and clears the flag. An edit made between runs becomes the new
 * baseline.
 *
 * ### Who takes part: excluded statuses, closing
 *
 * An agent in a status or stage of contact weight 0 (`set_contact()`) is out
 * of every layer. `exclude_status()` takes the agents in a status out of this
 * layer only: e.g., the school layer excludes "Quarantined", while the
 * household layer keeps them. An agent out of a layer neither samples it nor is
 * sampled through it, in both directions, and a mixing layer counts only the
 * agents present (its contact rate stays constant; within groups, the
 * probability per pair does). Who is out is the status, so the status history
 * already counts them. Excluded statuses are a setting, made between runs.
 *
 * `set_schedule()` makes a layer sample contacts on some days only, routinely
 * (work and transit on weekdays): a day off the schedule is a day closed, but
 * the schedule is a setting, made between runs, while `close()` is the
 * one-off intervention a global event makes during a run. `is_closed()` says
 * whether the layer has contacts today, for either reason.
 *
 * `close()` takes the whole layer out (a school closure) until `reopen()`:
 * it samples no contacts, so it logs none and transmits nothing, and the
 * days it was closed are not days shared in `layer::Mixing::groupmates_of()`.
 * Like any edit, closing during a run takes effect at the end of the step
 * (the layer is closed from the next day on) and lasts that run only; closing
 * between runs makes the layer closed in every run. A weight of 0
 * (`set_weight()`) is different: the contacts still happen (and are logged),
 * but transmit nothing.
 *
 * ### Contact tracing
 *
 * `enable_tracing()` makes the layer keep a log of its contacts (a
 * `ContactTracing`): every contact the layer samples, with a carrier, is
 * recorded under the carrier, with the day. The log starts empty at every run
 * and keeps the last `max_contacts` contacts of each agent. Tracing is off by
 * default, and it is a setting, not an edit: it lasts across runs whenever it
 * is turned on or off.
 *
 * @ingroup contact_layers
 */
class ContactLayer {

    friend class Model;

private:

    std::string name;
    size_t id = 0u;          ///< Position in the model's list, set by the model
    bool dirty = false;  ///< Edited since the baseline was set
    bool in_run = false; ///< Set by the model while it runs

    std::optional< ParamRef > weight_ref; ///< Empty: weight 1
    DiseaseTargets targets_;              ///< Diseases that spread here (see add_target())
    size_t trace_max = 0u;   ///< Contacts kept per agent; 0: not traced
    ContactTracing trace_log; ///< This run's contacts, if traced

    /// An empty log for a run of `n` agents (the model calls it at reset).
    void reset_tracing(size_t n)
    {
        if (trace_max > 0u)
            trace_log.reset(n, trace_max);
    }

    std::vector< std::string > excluded_statuses; ///< See exclude_status()
    std::vector< char > state_out;      ///< [state] Out of this layer by its statuses; empty: none
    std::vector< bool > state_silent;   ///< [state] Transmits nothing here; empty: the model's

    bool closed_base = false; ///< Closed in every run (closed between runs)
    bool closed = false;      ///< No contacts today (closed, or off its schedule)
    bool closed_next = false; ///< Closed from tomorrow (applied at the step's end)
    /// The spans of days without contacts this run, `[first, last]`
    /// (`INT_MAX`: still)
    std::vector< std::pair< int, int > > closed_spans;

    /// The days the layer samples contacts (empty: every day; see set_schedule())
    std::function< bool(int) > schedule;

    /// The masks of the run and the baseline closure (the model calls it at reset).
    void reset_participation(const Model & model);

    /// Applies today's `close()` / `reopen()` and the schedule, for tomorrow
    /// (the model calls it at the step's end).
    void apply_closure(int today);

protected:

    /// The model the layer belongs to (set by the model), or `nullptr`.
    Model * model = nullptr;

    /**
     * @brief Call it from every method that edits the layer.
     * @return `true` if the edit is to the baseline (made between runs);
     * otherwise it marks the layer dirty and returns `false`.
     */
    bool edits_baseline()
    {
        dirty = dirty || in_run;
        return !in_run;
    }

    /// Call it from `reset()` once the baseline has been restored.
    void mark_clean() { dirty = false; }

    /// Whether the model runs (settings are made between runs).
    bool in_a_run() const { return in_run; }

public:

    explicit ContactLayer(std::string name_) : name(std::move(name_)) {}

    /// A copy belongs to no model (see `clone()`).
    ContactLayer(const ContactLayer & other) :
        name(other.name), dirty(other.dirty), weight_ref(other.weight_ref),
        targets_(other.targets_), trace_max(other.trace_max),
        trace_log(other.trace_log),
        excluded_statuses(other.excluded_statuses), state_out(other.state_out),
        state_silent(other.state_silent), closed_base(other.closed_base),
        closed(other.closed), closed_next(other.closed_next),
        closed_spans(other.closed_spans), schedule(other.schedule) {}

    ContactLayer & operator=(const ContactLayer &) = delete;

    virtual ~ContactLayer() = default;

    const std::string & get_name() const { return name; }
    size_t get_id() const { return id; }
    bool is_dirty() const { return dirty; }

    /// The weight is the model parameter `param` (see the class description).
    void set_weight(std::string param) { weight_ref.emplace(std::move(param)); }

    bool has_weight() const { return weight_ref.has_value(); }

    /// Current weight in `model`.
    epiworld_double get_weight(const Model & model) const
    {
        return weight_ref ? (*weight_ref)(model) : 1.0;
    }

    /**
     * @name Disease targets
     * @brief Restrict the layer to the listed diseases, by name (see
     * `Tool::add_target()`, which these mirror). Without targets, every
     * disease spreads through the layer.
     */
    ///@{
    void add_target(std::string disease) { targets_.add(std::move(disease)); }
    const std::vector< std::string > & get_targets() const { return targets_.names; }
    ///@}

    /**
     * @name Contact tracing
     * @brief Keep a log of the layer's contacts (see the class description).
     */
    ///@{
    /// Turn tracing on, keeping up to `max_contacts` contacts per agent.
    void enable_tracing(size_t max_contacts = EPI_MAX_TRACKING);
    void disable_tracing();
    bool is_tracing() const { return trace_max > 0u; }
    /// The log of the current run (throws if tracing is off).
    ContactTracing & get_tracing();

    /// Log that `source` (the carrier) met `contact` today, if traced.
    void trace(size_t source, size_t contact)
    {
        if (trace_max > 0u)
            trace_log.add_contact(
                source, contact, static_cast< size_t >(model->today())
            );
    }

    /// Log that each of the `n` `sources` met `contact` today, if traced.
    void trace(const size_t * sources, size_t n, size_t contact)
    {
        if (trace_max > 0u)
            for (size_t k = 0u; k < n; ++k)
                trace(sources[k], contact);
    }
    ///@}

    /**
     * @name Who takes part
     * @brief Statuses out of this layer, and closing it (see the class
     * description).
     */
    ///@{
    /// Agents in `status` (a name in the model's statuses) are out of the layer.
    void exclude_status(std::string status);
    const std::vector< std::string > & get_excluded_statuses() const
    {
        return excluded_statuses;
    }

    void close();  ///< Closed from the next day (or every run, between runs)
    void reopen(); ///< Open from the next day (or every run, between runs)
    /// No contacts today: closed, or off its schedule. After a run, "today"
    /// is the day after its last (the run's last step prepares it).
    bool is_closed() const { return closed; }

    /**
     * @brief The days the layer samples contacts, as a cycle (between runs).
     * @param days Day `d` (contacts start on day 1) is on if `days[(d - 1) %
     * days.size()]`; e.g., five `true` and two `false` for weekdays only, the
     * run starting on a Monday.
     */
    void set_schedule(std::vector< bool > days);

    /// The days the layer samples contacts: those where `on(day)` is true
    /// (between runs; read at the end of the day before).
    void set_schedule(std::function< bool(int day) > on);

    void clear_schedule(); ///< Every day (the default; between runs)
    bool has_schedule() const { return static_cast< bool >(schedule); }

    /// Whether the schedule has day `day` on (always, without one); a closure
    /// takes the layer out on top of it.
    bool is_scheduled(int day) const { return !schedule || schedule(day); }

    /// Whether an agent in engine state `state` is out of the layer.
    bool excludes(epiworld_fast_uint state) const;

    /// Whether agent `agent_id` is out of this layer (by its excluded
    /// statuses, or a contact weight of 0), for samplers of their own.
    bool excludes_agent(size_t agent_id) const;

    /// The states out of this layer by its excluded statuses (`nullptr`: none;
    /// `excludes()` adds the states out of every layer).
    const std::vector< char > * own_exclusions() const
    {
        return state_out.empty() ? nullptr : &state_out;
    }

    /// Calls `f(first, last)` for each span of days in `[a, b]` the layer was open.
    template< typename F >
    void open_days(int a, int b, F f) const;
    ///@}

    /// Whether the virus spreads through the layer.
    bool admits(const Virus & v) const
    {
        return targets_.has(v);
    }

    /// Restore the baseline if dirty and size the buffers for a new run.
    virtual void reset(Model & model) = 0;

    /**
     * @brief Appends to `out` the contacts of agent `agent` this step that
     * carry a virus (repeats allowed), and reports them to the model's
     * post-sampling callback, if there is one.
     * @details This is the layer's part of a pull. Contacts in a state flagged
     * in `exclude` (if not null) are skipped: the model passes the layer's
     * `silent()` mask (states that transmit nothing here). Whether the contacts are able to
     * infect (their viruses, the weight) is up to the model.
     * @return The number of contacts appended.
     */
    virtual size_t sample(
        size_t agent,
        Model & model,
        const std::vector< bool > * exclude,
        std::vector< size_t > & out
    ) = 0;

    /// States that transmit nothing through the layer (not infectious, or
    /// out of it), or `nullptr` if every state transmits.
    const std::vector< bool > * silent() const;

    /// Rebuild what depends on the state of the population (end of each step).
    virtual void update(Model & model) = 0;

    /// Deep copy (the model's copy constructor calls it).
    virtual std::shared_ptr< ContactLayer > clone() const = 0;

    /// The layer as a network layer, or `nullptr` if it is of another kind.
    virtual layer::Network * as_network() { return nullptr; }

    /// The layer as a mixing layer, or `nullptr` if it is of another kind.
    virtual layer::Mixing * as_mixing() { return nullptr; }

    /**
     * @brief Whether the model's queue sees this layer's contacts.
     * @details True when a carrier can reach only agents the layer can list
     * from the carrier's side: its ties in an undirected network, its group in
     * a within-group mixing layer. A model with any other layer visits every
     * agent instead of using the queue.
     */
    virtual bool supports_queue() const { return false; }

};

inline void ContactLayer::enable_tracing(size_t max_contacts)
{
    if (max_contacts == 0u)
        throw std::invalid_argument(
            "Contact tracing should keep at least one contact per agent."
        );

    trace_max = max_contacts;

    // Turned on in a model that has agents: logs from now on
    if (model != nullptr)
        reset_tracing(model->size());
}

inline void ContactLayer::disable_tracing()
{
    trace_max = 0u;
    trace_log = ContactTracing();
}

inline ContactTracing & ContactLayer::get_tracing()
{
    if (trace_max == 0u)
        throw std::logic_error(
            "The contact layer '" + name + "' is not traced; call "
            "enable_tracing() first."
        );

    return trace_log;
}

#endif
