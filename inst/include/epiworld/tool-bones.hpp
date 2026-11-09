#ifndef EPIWORLD_TOOL_BONES_HPP
#define EPIWORLD_TOOL_BONES_HPP

class Virus;

class Agent;

class Model;

class Tool;

/**
 * @brief Tools for defending the agent against the virus
 * 
 */
class Tool {
    friend class Agent;
    friend class Model;
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


    /// Bitmask of targeted virus lineages (all ones: every virus).
    uint64_t target_mask = ~uint64_t(0);


    void set_agent(Agent * p, size_t idx);

public:

    /// Bit of a virus lineage in a target mask (contact layers use it too).
    static uint64_t target_bit(int lineage_id); ///< Throws for ids out of range.

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
     * @name Virus targets
     * 
     * @details
     * By default, a tool acts on every virus. Adding targets restricts the
     * tool (all four effects: susceptibility and transmission reduction,
     * recovery enhancer, and death reduction) to the listed virus lineages.
     * A lineage is a virus added with `Model::add_virus()` together with all
     * of its mutations, so a tool targeting a virus also acts on its
     * variants. Lineage ids are the virus ids assigned by
     * `Model::add_virus()`; only lineages 0 to 62 can be targeted.
     * 
     * @param lineage_id Id of the virus lineage (see
     * `Virus::get_lineage_id()`).
     * @param v A virus already added to the model.
     */
    ///@{
    void add_target(int lineage_id);
    void add_target(const Virus & v);
    void set_targets(const std::vector< int > & lineage_ids);
    std::vector< int > get_targets() const; ///< Empty if the tool acts on every virus.
    void clear_targets();
    bool targets(const Virus & v) const;
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