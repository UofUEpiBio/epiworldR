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
        name(other.name), dirty(other.dirty) {}

    ContactLayer<TSeq> & operator=(const ContactLayer<TSeq> &) = delete;

    virtual ~ContactLayer() = default;

    const std::string & get_name() const { return name; }
    size_t get_id() const { return id; }
    bool is_dirty() const { return dirty; }

    /// Restore the baseline if dirty and size the buffers for a new run.
    virtual void reset(Model<TSeq> & model) = 0;

    /// Rebuild what depends on the state of the population (end of each step).
    virtual void update(Model<TSeq> & model) = 0;

    /// Deep copy (the model's copy constructor calls it).
    virtual std::shared_ptr< ContactLayer<TSeq> > clone() const = 0;

    /// The layer as a network layer, or `nullptr` if it is of another kind.
    virtual layer::Network<TSeq> * as_network() { return nullptr; }

};

#endif
