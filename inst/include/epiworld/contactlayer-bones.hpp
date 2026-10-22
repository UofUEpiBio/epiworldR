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
 *   (entities, viruses, tools, the global events' set-up), so that those see
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
 * off). `add_virus_target()` restricts the layer to some viruses, as tools are
 * restricted; by default every virus spreads through every layer.
 *
 * ### Edits and the dirty flag
 *
 * Whether an edit is permanent depends only on when it is made. An edit made
 * during a run (from an update function or a global event) applies to that run
 * only: it sets the *dirty* flag, and `reset()` restores the baseline of a
 * dirty layer and clears the flag. An edit made between runs becomes the new
 * baseline.
 *
 * @ingroup contact_layers
 */
template<typename TSeq>
class ContactLayer {

    friend class Model<TSeq>;

private:

    std::string name;
    size_t id = 0u;          ///< Position in the model's list, set by the model
    bool dirty = false;  ///< Edited since the baseline was set
    bool in_run = false; ///< Set by the model while it runs

    std::optional< ParamRef > weight_ref; ///< Empty: weight 1
    uint64_t virus_mask = ~uint64_t(0);   ///< Lineages that spread here

protected:

    /// The model the layer belongs to (set by the model), or `nullptr`.
    Model<TSeq> * model = nullptr;

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

public:

    explicit ContactLayer(std::string name_) : name(std::move(name_)) {}

    /// A copy belongs to no model (see `clone()`).
    ContactLayer(const ContactLayer<TSeq> & other) :
        name(other.name), dirty(other.dirty), weight_ref(other.weight_ref),
        virus_mask(other.virus_mask) {}

    ContactLayer<TSeq> & operator=(const ContactLayer<TSeq> &) = delete;

    virtual ~ContactLayer() = default;

    const std::string & get_name() const { return name; }
    size_t get_id() const { return id; }
    bool is_dirty() const { return dirty; }

    /// The weight is the model parameter `param` (see the class description).
    void set_weight(std::string param) { weight_ref.emplace(std::move(param)); }

    bool has_weight() const { return weight_ref.has_value(); }

    /// Current weight in `model`.
    epiworld_double get_weight(const Model<TSeq> & model) const
    {
        return weight_ref ? (*weight_ref)(model) : 1.0;
    }

    /**
     * @name Virus targets
     * @brief Restrict the layer to the listed virus lineages (see
     * `Tool::add_target()`, which these mirror). The first target turns the
     * default (every virus) off.
     */
    ///@{
    void add_virus_target(int lineage_id);
    void add_virus_target(const Virus<TSeq> & v);
    ///@}

    /// Whether the virus spreads through the layer.
    bool admits(const Virus<TSeq> & v) const
    {
        return (virus_mask & v.lineage_bit) != 0u;
    }

    /// Restore the baseline if dirty and size the buffers for a new run.
    virtual void reset(Model<TSeq> & model) = 0;

    /**
     * @brief Appends to `out` the contacts of agent `agent` this step that
     * carry a virus (repeats allowed), and reports them to the model's
     * post-sampling callback, if there is one.
     * @details This is the layer's part of a pull. Contacts in a state flagged
     * in `exclude` (if not null) are skipped. Whether the contacts are able to
     * infect (their viruses, the weight) is up to the model.
     * @return The number of contacts appended.
     */
    virtual size_t sample(
        size_t agent,
        Model<TSeq> & model,
        const std::vector< bool > * exclude,
        std::vector< size_t > & out
    ) = 0;

    /// Rebuild what depends on the state of the population (end of each step).
    virtual void update(Model<TSeq> & model) = 0;

    /// Deep copy (the model's copy constructor calls it).
    virtual std::shared_ptr< ContactLayer<TSeq> > clone() const = 0;

    /// The layer as a network layer, or `nullptr` if it is of another kind.
    virtual layer::Network<TSeq> * as_network() { return nullptr; }

};

template<typename TSeq>
inline void ContactLayer<TSeq>::add_virus_target(int lineage_id)
{

    const uint64_t bit = Tool<TSeq>::target_bit(lineage_id);

    // The first target turns the default off
    if (virus_mask == ~uint64_t(0))
        virus_mask = 0u;

    virus_mask |= bit;

}

template<typename TSeq>
inline void ContactLayer<TSeq>::add_virus_target(const Virus<TSeq> & v)
{

    if (v.get_lineage_id() < 0)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no lineage id. Add it to "
            "the model with Model::add_virus() before using it as a target."
        );

    add_virus_target(v.get_lineage_id());

}

#endif
