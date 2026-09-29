#ifndef EPIWORLD_PARAM_REF_HPP
#define EPIWORLD_PARAM_REF_HPP

/**
 * @brief Position of a parameter in a model's parameter table.
 *
 * @details Returned by `Model::get_param_id()` and read with
 * `Model::par_at()`. Parameters are never removed, so an id stays valid for
 * the model that returned it and for every copy of that model (including the
 * copies `run_multiple()` makes). Models built separately may place the same
 * name at different positions; use `ParamRef` when the same code runs on
 * different models.
 */
struct ParamId {
    size_t idx;
};

/**
 * @brief Identifies one layout of a model's parameter table.
 *
 * @details Every model starts with a fresh id, copies keep the id of the
 * model they copy, and adding a new parameter gives the model a fresh id. So
 * two models with the same id always map names to the same positions. Zero
 * is never returned; `ParamRef` uses it for "not resolved yet".
 */
inline uint32_t new_param_layout_id()
{
    static std::atomic< uint32_t > counter{0u};
    uint32_t id = ++counter;
    while (id == 0u)
        id = ++counter;
    return id;
}

/**
 * @brief A model parameter referenced by name, resolved to its position on
 * first use.
 *
 * @details Calling the object with a model returns the parameter's current
 * value. The first call on a model looks the name up in the parameter table
 * and caches its position together with the model's parameter layout (see
 * `new_param_layout_id()`); later calls on that model, or on any copy of it,
 * read the value directly. Using the object with a model of a different
 * layout looks the name up again, so a `ParamRef` is always safe to share
 * between models and threads.
 *
 * Values are never cached, only positions: changes made with `set_param()`
 * are seen immediately.
 *
 * `EPI_PAR(model, "name")` wraps a function-local `static ParamRef`, which is
 * the easiest way to use it inside update functions.
 */
class ParamRef {
private:
    std::string pname;

    // (layout id << 32) | position; 0 means not resolved.
    mutable std::atomic< uint64_t > cache{0u};

public:

    explicit ParamRef(std::string name) : pname(std::move(name)) {};

    ParamRef(const ParamRef & other) :
        pname(other.pname),
        cache(other.cache.load(std::memory_order_relaxed)) {};

    ParamRef & operator=(const ParamRef & other)
    {
        pname = other.pname;
        cache.store(
            other.cache.load(std::memory_order_relaxed),
            std::memory_order_relaxed
        );
        return *this;
    };

    const std::string & name() const { return pname; };

    /**
     * @brief Position of the parameter in `model`'s table.
     * @throws std::logic_error if the model has no parameter with this name.
     */
    template< typename TModel >
    ParamId id(const TModel & model) const
    {
        const uint32_t layout = model.get_param_layout_id();
        const uint64_t c = cache.load(std::memory_order_relaxed);
        if (static_cast< uint32_t >(c >> 32) == layout)
            return ParamId{static_cast< size_t >(c & 0xffffffffu)};

        ParamId res = model.get_param_id(pname);
        if (res.idx > 0xffffffffu)
            return res;

        cache.store(
            (static_cast< uint64_t >(layout) << 32) |
                static_cast< uint64_t >(res.idx),
            std::memory_order_relaxed
        );
        return res;
    };

    /**
     * @brief Current value of the parameter in `model`.
     * @throws std::logic_error if the model has no parameter with this name.
     */
    template< typename TModel >
    epiworld_double operator()(const TModel & model) const
    {
        return model.par_at(id(model));
    };

};

#endif
