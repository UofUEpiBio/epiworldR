#ifndef EPIWORLD_TOOL_BONES_HPP
#define EPIWORLD_TOOL_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;

template<typename TSeq>
class Model;

template<typename TSeq>
class Tool;

/**
 * @brief Tools for defending the agent against the virus
 * 
 * @tparam TSeq Type of sequence
 */
template<typename TSeq> 
class Tool {
    friend class Agent<TSeq>;
    friend class Model<TSeq>;
protected:

    Agent<TSeq> * agent = nullptr;
    int pos_in_agent        = -99; ///< Location in the agent

    int date = -99;
    int id   = -99;
    std::string tool_name;
    
    EPI_TYPENAME_TRAITS(TSeq, int) sequence = 
        EPI_TYPENAME_TRAITS(TSeq, int)(); ///< Sequence of the tool

    ToolFun<TSeq> susceptibility_reduction = nullptr;
    ToolFun<TSeq> transmission_reduction   = nullptr;
    ToolFun<TSeq> recovery_enhancer        = nullptr;
    ToolFun<TSeq> death_reduction          = nullptr;

    ToolToAgentFun<TSeq> dist = nullptr;

    /// Bitmask of targeted virus lineages (all ones: every virus).
    uint64_t target_mask = ~uint64_t(0);
    static uint64_t target_bit(int lineage_id); ///< Throws for ids out of range.

    epiworld_fast_int state_init = -99;
    epiworld_fast_int state_post = -99;

    epiworld_fast_int queue_init = Queue<TSeq>::NoOne; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = Queue<TSeq>::NoOne; ///< Change of state when removed from agent.

    void set_agent(Agent<TSeq> * p, size_t idx);

public:
    Tool();
    Tool(std::string name = "unknown tool");
    Tool(
        std::string name,
        epiworld_double prevalence,
        bool as_proportion
    );

    virtual ~Tool() = default;

    void set_sequence(TSeq d);
    void set_sequence(std::shared_ptr<TSeq> d);
    EPI_TYPENAME_TRAITS(TSeq, int) get_sequence();

    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    virtual epiworld_double get_susceptibility_reduction(VirusPtr<TSeq> & v, Model<TSeq> * model);
    virtual epiworld_double get_transmission_reduction(VirusPtr<TSeq> & v, Model<TSeq> * model);
    virtual epiworld_double get_recovery_enhancer(VirusPtr<TSeq> & v, Model<TSeq> * model);
    virtual epiworld_double get_death_reduction(VirusPtr<TSeq> & v, Model<TSeq> * model);
    
    virtual void set_susceptibility_reduction_fun(ToolFun<TSeq> fun);
    virtual void set_transmission_reduction_fun(ToolFun<TSeq> fun);
    virtual void set_recovery_enhancer_fun(ToolFun<TSeq> fun);
    virtual void set_death_reduction_fun(ToolFun<TSeq> fun);

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
    void add_target(const Virus<TSeq> & v);
    void set_targets(const std::vector< int > & lineage_ids);
    std::vector< int > get_targets() const; ///< Empty if the tool acts on every virus.
    void clear_targets();
    bool targets(const Virus<TSeq> & v) const;
    ///@}

    void set_name(std::string name);
    virtual std::string get_name() const;

    Agent<TSeq> * get_agent();
    int get_id() const;
    void set_id(int id);
    void set_date(int d);
    int get_date() const;

    void set_state(epiworld_fast_int init, epiworld_fast_int post);
    void set_queue(epiworld_fast_int init, epiworld_fast_int post);
    void get_state(epiworld_fast_int * init, epiworld_fast_int * post);
    void get_queue(epiworld_fast_int * init, epiworld_fast_int * post);

    bool operator==(const Tool<TSeq> & other) const;
    bool operator!=(const Tool<TSeq> & other) const {return !operator==(other);};

    void print() const;

    void distribute(Model<TSeq> * model);
    void set_distribution(ToolToAgentFun<TSeq> fun);

    virtual std::unique_ptr<Tool<TSeq>> clone_ptr() const; 

};

#endif