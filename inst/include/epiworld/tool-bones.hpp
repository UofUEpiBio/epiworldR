#ifndef EPIWORLD_TOOL_BONES_HPP
#define EPIWORLD_TOOL_BONES_HPP

class Virus;

class Agent;

class Model;

class Tool;

/**
 * @brief The diseases a tool or a contact layer acts on, by name.
 * @details Resolved to a mask of the model's diseases when a run starts
 * (`resolve()`, see `Model::reset()`); all of them until then.
 */
struct DiseaseTargets {
    std::vector< std::string > names; ///< Empty: see `resolve()`
    uint64_t mask = ~uint64_t(0);     ///< Bit `d`: disease `d`

    void add(std::string disease);

    /**
     * Makes the names a mask. Without names: every disease, unless `required`
     * and the model has several (an error naming `owner`).
     */
    void resolve(const Model & model, const std::string & owner, bool required);

    /// Whether it acts on `v` (a virus not in a model: no)
    bool has(const Virus & v) const;

    bool operator==(const DiseaseTargets & o) const { return names == o.names; }
};

/**
 * @brief The contact layers a tool's contact effects act in, by name (see
 * `Tool::add_layer_target()`).
 */
struct LayerTargets {
    std::vector< std::string > names; ///< Empty: every layer
    uint64_t mask = ~uint64_t(0);     ///< Bit `l`: layer `l` (once resolved)
    bool resolved = false;            ///< Whether `mask` matches `names`

    void add(std::string layer);

    /// Makes the names a mask (an unknown layer is an error naming `owner`).
    void resolve(Model & model, const std::string & owner);

    /// Whether it acts through layer `layer` (-1: not through a layer, where
    /// only a tool without layer targets acts). With names, `resolve()` must
    /// have run since they last changed (`Tool::acts_in()` checks).
    bool has(int layer) const;

    bool operator==(const LayerTargets & o) const { return names == o.names; }
};

/**
 * @brief Tools for defending the agent against the virus
 * 
 */
class Tool {
    friend class Agent;
    friend class Model;
    friend class DataBase;
protected:

    Agent * agent = nullptr;
    int pos_in_agent        = -99; ///< Location in the agent

    int date = -99;
    int id   = -99;
    std::string tool_name;
    

    ToolFun susceptibility_reduction = nullptr;
    ToolFun transmission_reduction   = nullptr;
    ToolFun recovery_enhancer        = nullptr;
    ToolFun death_reduction          = nullptr;


    DiseaseTargets targets_; ///< The diseases it acts on (see add_target())
    LayerTargets layers_;    ///< The layers its contact effects act in


    void set_agent(Agent * p, size_t idx);

public:

    Tool();
    Tool(std::string name = "unknown tool");

    virtual ~Tool() = default;


    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    virtual epiworld_double get_susceptibility_reduction(VirusPtr & v, Model * model);
    virtual epiworld_double get_transmission_reduction(VirusPtr & v, Model * model);
    virtual epiworld_double get_recovery_enhancer(VirusPtr & v, Model * model);
    virtual epiworld_double get_death_reduction(VirusPtr & v, Model * model);
    
    virtual void set_susceptibility_reduction_fun(ToolFun fun);
    virtual void set_transmission_reduction_fun(ToolFun fun);
    virtual void set_recovery_enhancer_fun(ToolFun fun);
    virtual void set_death_reduction_fun(ToolFun fun);

    virtual void set_susceptibility_reduction(std::string param);
    virtual void set_transmission_reduction(std::string param);
    virtual void set_recovery_enhancer(std::string param);
    virtual void set_death_reduction(std::string param);

    // Deleting pointer versions to avoid mistakes
    virtual void set_susceptibility_reduction(epiworld_double * prob) = delete;
    virtual void set_transmission_reduction(epiworld_double * prob) = delete;
    virtual void set_recovery_enhancer(epiworld_double * prob) = delete;
    virtual void set_death_reduction(epiworld_double * prob) = delete;

    virtual void set_susceptibility_reduction(epiworld_double prob);
    virtual void set_transmission_reduction(epiworld_double prob);
    virtual void set_recovery_enhancer(epiworld_double prob);
    virtual void set_death_reduction(epiworld_double prob);
    ///@}

    /**
     * @name Disease targets
     *
     * @details A tool acts (all four effects: susceptibility and transmission
     * reduction, recovery enhancer, and death reduction) on the viruses of
     * the diseases it targets, every strain and variant of each, by the
     * disease's name (see `Model::get_diseases()`). In a model with one
     * disease, a tool without targets acts on it; with several, a tool must
     * name its diseases (an error at the start of the run otherwise): a
     * measles vaccine does not protect against influenza unless it says so.
     *
     * @param disease The name of a disease of the model.
     */
    ///@{
    void add_target(std::string disease);
    void set_targets(std::vector< std::string > diseases);
    const std::vector< std::string > & get_targets() const { return targets_.names; }
    void clear_targets() { targets_ = DiseaseTargets(); }
    bool targets(const Virus & v) const { return targets_.has(v); }
    ///@}

    /**
     * @name Layer targets
     *
     * @details A tool's contact effects (susceptibility and transmission
     * reduction) act on the contacts through the layers it targets, by the
     * layer's name (see `Model::add_layer()`): e.g., a mask worn at work and
     * school, not at home. Without layer targets, they act in every layer.
     * The recovery enhancer and death reduction are not about contacts, so
     * layer targets leave them alone. An infection not drawn from a layer
     * (e.g., given by code) sees only the tools without layer targets. A
     * name the model has no layer of is an error at the start of the run; at
     * most 64 layers can be targeted (ids 0 to 63).
     *
     * @param layer The name of a layer of the model.
     */
    ///@{
    void add_layer_target(std::string layer);
    void set_layer_targets(std::vector< std::string > layers);
    const std::vector< std::string > & get_layer_targets() const { return layers_.names; }
    void clear_layer_targets() { layers_ = LayerTargets(); }
    /// Whether its contact effects act through layer `layer` (-1: no layer).
    /// With layer targets, an error until they are resolved: when the run
    /// starts, or when the tool is added to a running model.
    bool acts_in(int layer) const;
    ///@}

    void set_name(std::string name);
    virtual std::string get_name() const;

    Agent * get_agent();
    int get_id() const;
    void set_id(int id);
    void set_date(int d);
    int get_date() const;


    bool operator==(const Tool & other) const;
    bool operator!=(const Tool & other) const {return !operator==(other);};

    void print() const;


    virtual std::unique_ptr<Tool> clone_ptr() const; 

};

#endif