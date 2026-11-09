#ifndef EPIWORLD_REGISTRY_BONES_HPP
#define EPIWORLD_REGISTRY_BONES_HPP

/**
 * @brief One argument of a registered model's constructor.
 * @ingroup epi_models
 *
 * @details `arg` is the argument's name, the one every wrapper uses (R, Python
 * and JS take `transmission_prob`, not each their own). `name` is the model
 * parameter the argument sets, as `Model::get_param()` knows it
 * ("Transmission probability"); it is empty for the arguments that are not
 * parameters (the population size `n`, the prevalence). `min` and `max` are
 * hard limits, checked by `ModelRegistry::make_model()` before the model is
 * built.
 *
 * The static functions build the kinds that repeat across models.
 */
struct ParamInfo {
    std::string arg;           ///< Argument name, e.g. "transmission_prob".
    std::string name;          ///< Model parameter it sets, or empty.
    double value;              ///< Default value.
    double min;                ///< Smallest valid value.
    double max;                ///< Largest valid value (may be infinite).
    bool integer;              ///< Whether the value must be a whole number.
    std::string description;   ///< One sentence, for help pages and UIs.

    /// A probability or proportion, in [0, 1].
    static ParamInfo probability(
        std::string arg, std::string name, double value,
        std::string description
    );

    /// A mean duration in days, in [0, inf).
    static ParamInfo duration(
        std::string arg, std::string name, double value,
        std::string description
    );

    /// A whole number of days, where a negative value (-1) turns it off.
    static ParamInfo period(
        std::string arg, std::string name, double value,
        std::string description
    );

    /// A whole number of agents or days, in [0, max].
    static ParamInfo count(
        std::string arg, std::string name, double value,
        std::string description
    );

    /// The population size `n`, at least 1.
    static ParamInfo agents(double value);

    /// The proportion of agents infected at the start, `prevalence`.
    static ParamInfo prevalence(double value);
};

/**
 * @brief The values a registered model's builder receives.
 * @ingroup epi_models
 *
 * @details Complete (every argument of the model, defaults filled in) and
 * checked against the model's `ParamInfo`. `a("transmission_prob")` reads a
 * value; `a.virus_name()` the name of the model's virus.
 */
class ModelArgs {
private:
    std::map< std::string, double > values;
    std::string vname;

public:
    ModelArgs(std::map< std::string, double > values_, std::string vname_);

    /// @throws std::logic_error if `arg` is not one of the model's arguments.
    double operator()(const std::string & arg) const;

    const std::string & virus_name() const noexcept;
};

/**
 * @brief A model the registry can build: its arguments and builder.
 * @ingroup epi_models
 *
 * @details Each registered model returns its own from a static
 * `model_info()`, next to its constructor. `params` lists the constructor's
 * arguments in its order. The builder calls the constructor and nothing else:
 * the model it returns has no layers, which are added the same way as for a
 * model built directly (`Model::add_layer()`). A virus a user adds to the
 * model is a strain of the model's own: a copy of `get_virus(0)` with its own
 * name and, if need be, its own parameters (`Virus::use_param()`), so it has
 * the same stages.
 */
struct ModelInfo {
    std::string id;                   ///< e.g. "SEIR", for make_model()
    std::string label;                ///< e.g. "Susceptible-Exposed-Infected-Recovered"
    std::string description;          ///< One paragraph
    std::string virus_name;           ///< Default virus name; empty if the model names its own
    std::vector< ParamInfo > params;  ///< Constructor arguments, in order
    std::function< std::unique_ptr< Model >(const ModelArgs &) > build;
};

/**
 * @brief Builds models by name, and lists them with their arguments.
 * @ingroup epi_models
 *
 * @details One place per model: the wrappers bind `list_models()` and
 * `make_model()` once instead of every constructor, and read argument names,
 * defaults, limits and descriptions from `ParamInfo`.
 *
 * A registry starts empty. Each library adds its models with its own
 * `register_models()`, defined in its umbrella header next to the includes of
 * the models it lists; a model's entry is its static `model_info()`.
 *
 * ```cpp
 * ModelRegistry registry;
 * epimodels::register_models(registry);            // the core models
 * measles::register_models(registry);              // and, optionally, measles'
 *
 * auto model = registry.make_model(
 *     "SEIR", {{"n", 5000}, {"transmission_prob", 0.1}}, "flu"
 * );                                               // the rest take defaults
 * model->add_layer(layer::network("contacts", graph::smallworld(10, 0.1)));
 * model->run(100, 123);
 * ```
 */
class ModelRegistry {
private:
    std::vector< ModelInfo > models;

public:
    /**
     * @brief Adds a model.
     * @throws std::invalid_argument if the id is empty or already taken, an
     * argument is repeated, a default is outside its limits, or there is no
     * builder.
     */
    ModelRegistry & add(ModelInfo info);

    /// Every model, in the order they were added.
    const std::vector< ModelInfo > & list_models() const noexcept;

    /// @throws std::invalid_argument naming the known ids.
    const ModelInfo & get_model_info(const std::string & id) const;

    /**
     * @brief Builds a model from its id and the values of its arguments.
     * @param id The model, as in `list_models()`.
     * @param values Argument values by `ParamInfo::arg`; the arguments left out
     * take their defaults.
     * @param virus_name Name of the model's virus; empty for the default. Only
     * for the models whose `ModelInfo::virus_name` is not empty.
     * @throws std::invalid_argument for an unknown id or argument, a value
     * outside its limits or not whole when it must be, or a virus name for a
     * model that names its own virus.
     */
    std::unique_ptr< Model > make_model(
        const std::string & id,
        const std::map< std::string, double > & values = {},
        const std::string & virus_name = ""
    ) const;
};

#endif
