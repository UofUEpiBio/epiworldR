#ifndef EPIWORLD_PERSON_BONES_HPP
#define EPIWORLD_PERSON_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Viruses;

template<typename TSeq>
class Viruses_const;

template<typename TSeq>
class Tool;

template<typename TSeq>
class Tools;

template<typename TSeq>
class Tools_const;

template<typename TSeq>
class Queue;

template<typename TSeq>
struct Event;

template<typename TSeq>
class Entity;

template<typename TSeq>
class Entities;

template<typename TSeq>
inline void default_add_virus(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_tool(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_entity(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_virus(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_tool(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_entity(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_change_state(Event<TSeq> & a, Model<TSeq> * m);



/**
 * @brief Agent (agents)
 * 
 * @tparam TSeq Sequence type (should match `TSeq` across the model)
 */
template<typename TSeq>
class Agent {
    friend class Model<TSeq>;
    friend class Virus<TSeq>;
    friend class Tool<TSeq>;
    friend class Tools<TSeq>;
    friend class Tools_const<TSeq>;
    friend class Queue<TSeq>;
    friend class Entities<TSeq>;
    friend class AgentsSample<TSeq>;
    friend void default_add_virus<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_add_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_add_entity<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_entity<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_change_state<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
private:
    
    Model<TSeq> * model;

    std::vector< size_t > * neighbors = nullptr;
    std::vector< size_t > * neighbors_locations = nullptr;
    size_t n_neighbors = 0u;

    std::vector< size_t > entities;
    std::vector< size_t > entities_locations;
    size_t n_entities = 0u;

    unsigned int state = 0u;
    unsigned int state_prev = 0u; ///< For accounting, if need to undo a change.
    
    int state_last_changed = -1; ///< Last time the agent was updated.
    int id = -1;
    
    VirusPtr<TSeq> virus = nullptr;

    std::vector< ToolPtr<TSeq> > tools;
    unsigned int n_tools = 0u;

public:

    Agent();
    Agent(Agent<TSeq> && p);
    Agent(const Agent<TSeq> & p);
    Agent<TSeq> & operator=(const Agent<TSeq> & other_agent);
    ~Agent();

    /**
     * @name Add/Remove Virus/Tool
     * 
     * Any of these is ultimately reflected at the end of the iteration.
     * 
     * @param tool Tool to add
     * @param virus Virus to add
     * @param state_new state after the change
     * @param queue 
     */
    ///@{
    void add_tool(
        ToolPtr<TSeq> tool,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_tool(
        Tool<TSeq> tool,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
        );

    void set_virus(
        VirusPtr<TSeq> virus,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
        );

    void set_virus(
        Virus<TSeq> virus,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
        );

    void add_entity(
        Entity<TSeq> & entity,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
        );

    void rm_tool(
        epiworld_fast_uint tool_idx,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_tool(
        ToolPtr<TSeq> & tool,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_virus(
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_entity(
        epiworld_fast_uint entity_idx,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_entity(
        Entity<TSeq> & entity,
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    );

    void rm_agent_by_virus(
        Model<TSeq> * model,
        epiworld_fast_int state_new = -99,
        epiworld_fast_int queue = -99
    ); ///< Agent removed by virus
    ///@}
    
    /**
     * @name Get the rates (multipliers) for the agent
     * 
     * @param v A pointer to a virus.
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_susceptibility_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_transmission_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_recovery_enhancer(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_death_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    ///@}

    int get_id() const; ///< Id of the individual

    VirusPtr<TSeq> & get_virus();
    const VirusPtr<TSeq> & get_virus() const;

    ToolPtr<TSeq> & get_tool(int i);
    Tools<TSeq> get_tools();
    const Tools_const<TSeq> get_tools() const;
    size_t get_n_tools() const noexcept;

    void mutate_virus();
    void add_neighbor(
        Agent<TSeq> & p,
        bool check_source = true,
        bool check_target = true
        );

    /**
     * @brief Swaps neighbors between the current agent and agent `other`
     * 
     * @param other 
     * @param n_this 
     * @param n_other 
     */
    void swap_neighbors(
        Agent<TSeq> & other,
        size_t n_this,
        size_t n_other
    );

    std::vector< Agent<TSeq> * > get_neighbors();
    size_t get_n_neighbors() const;

    void change_state(
        Model<TSeq> * model,
        epiworld_fast_uint new_state,
        epiworld_fast_int queue = 0
        );

    const unsigned int & get_state() const;

    void reset();

    bool has_tool(epiworld_fast_uint t) const;
    bool has_tool(std::string name) const;
    bool has_tool(const Tool<TSeq> & t) const;
    bool has_virus(epiworld_fast_uint t) const;
    bool has_virus(std::string name) const;
    bool has_virus(const Virus<TSeq> & v) const;
    bool has_entity(epiworld_fast_uint t) const;
    bool has_entity(std::string name) const;

    void print(Model<TSeq> * model, bool compressed = false) const;

    /**
     * @brief Access the j-th column of the agent
     * 
     * If an external array has been specified, then these two
     * functions can be used to access additional agent's features 
     * not included in the model.
     * 
     * The `operator[]` method is with no boundary check, whereas
     * the `operator()` method checks boundaries. The former can result
     * in a segfault.
     * 
     * 
     * @param j 
     * @return double& 
     */
    ///@{
    double & operator()(size_t j);
    double & operator[](size_t j);
    double operator()(size_t j) const;
    double operator[](size_t j) const;
    ///@}

    Entities<TSeq> get_entities();
    const Entities_const<TSeq> get_entities() const;
    const Entity<TSeq> & get_entity(size_t i) const;
    Entity<TSeq> & get_entity(size_t i);
    size_t get_n_entities() const;

    bool operator==(const Agent<TSeq> & other) const;
    bool operator!=(const Agent<TSeq> & other) const {return !operator==(other);};

};



#endif