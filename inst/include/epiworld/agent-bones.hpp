#ifndef EPIWORLD_PERSON_BONES_HPP
#define EPIWORLD_PERSON_BONES_HPP

class Model;

class Virus;

class Tool;


class Queue;

struct Event;

class AgentsSample;

/**
 * @brief Non-allocating range over an agent's neighbors in a network layer.
 *
 * @details It iterates the agent's neighbor ids in place and resolves each one
 * against the model's population as it goes, so a loop over the neighbors does
 * not allocate:
 *
 * ```cpp
 * for (auto * neighbor : net.neighbors_view(p->get_id()))
 *     ...
 * ```
 *
 * The order is the layer's tie order (see `layer::Network`). The view borrows
 * from the layer and the model, so it must not outlive either, and it is
 * invalidated by anything that changes the agent's ties.
 *
 */
class NeighborsView {
private:

    const size_t * first = nullptr;
    size_t n = 0u;
    std::vector< Agent > * pop = nullptr;

public:

    class iterator {

        friend class NeighborsView;

        const size_t * ptr = nullptr;
        std::vector< Agent > * pop = nullptr;

        iterator(const size_t * ptr, std::vector< Agent > * pop) :
            ptr(ptr), pop(pop) {}

    public:

        Agent * operator*() const { return &pop->operator[](*ptr); }
        iterator & operator++() { ++ptr; return *this; }
        bool operator!=(const iterator & other) const { return ptr != other.ptr; }
        bool operator==(const iterator & other) const { return ptr == other.ptr; }

    };

    NeighborsView() = default;

    NeighborsView(
        const size_t * first,
        size_t n,
        std::vector< Agent > * pop
    ) : first(first), n(n), pop(pop) {}

    // The size is carried rather than derived as `last - first`: an agent with
    // no ties has nothing to point at, and subtracting two null pointers is
    // undefined behaviour (they point into no array).
    iterator begin() const { return iterator(first, pop); }
    iterator end() const { return iterator(first + n, pop); }

    size_t size() const { return n; }
    bool empty() const { return n == 0u; }

    /// Id of the `k`-th neighbor, without touching the neighbor's `Agent`.
    /// @throws std::out_of_range if `k >= size()`.
    size_t id(size_t k) const
    {
        if (k >= n)
            throw std::out_of_range("NeighborsView::id(): index out of range.");
        return first[k];
    }

    /// The `k`-th neighbor.
    /// @throws std::out_of_range if `k >= size()`.
    Agent * agent(size_t k) const
    {
        if (k >= n)
            throw std::out_of_range("NeighborsView::agent(): index out of range.");
        return &pop->operator[](first[k]);
    }

    /// The ids of the `size()` neighbors, in order (like `std::vector::data()`:
    /// only valid for indices below `size()`; `nullptr` when there are none).
    /// For scans that must not touch the neighbors' `Agent`s.
    const size_t * ids() const { return first; }

};

/**
 * @brief Agent (agents)
 * 
 */
class Agent {
    friend class Model;
    friend class Virus;
    friend class Tool;
    friend class Queue;
    friend class AgentsSample;
    friend class DataBase;
    friend class layer::Network;
protected:

    unsigned int state = 0u;
    unsigned int state_prev = 0u; ///< For accounting, if need to undo a change.
    
    int state_last_changed = -1; ///< Last time the agent was updated.
    int id = -1;
    
    VirusPtr virus = nullptr;

    std::vector< ToolPtr > tools;

    void reset(); ///< Resets the agent to the initial state (no virus, no tools, state 0.)

    /**
     * @name The engine's changes (with the target state and queue)
     * @details -99 means the default: no change of state (tools), the
     * infection's target (a virus given from outside).
     */
    ///@{
    void add_tool(
        Model & model,
        const Tool & tool,
        epiworld_fast_int state_new,
        epiworld_fast_int queue
        );

    void set_virus(
        Model & model,
        const Virus & virus,
        epiworld_fast_int state_new,
        epiworld_fast_int queue
        );

    void rm_tool(
        Model & model,
        epiworld_fast_uint tool_idx,
        epiworld_fast_int state_new,
        epiworld_fast_int queue
    );

    void rm_tool(
        Model & model,
        ToolPtr & tool,
        epiworld_fast_int state_new,
        epiworld_fast_int queue
    );

    void rm_virus(
        Model & model,
        epiworld_fast_int state_new,
        epiworld_fast_int queue
    );

    void change_state(
        Model & model,
        epiworld_fast_uint new_state,
        epiworld_fast_int queue = 0
        );

    /// The engine's index of the agent's stage and status (see
    /// `Model::get_stage()` and `Model::get_status()`).
    unsigned int get_state() const;
    unsigned int get_state_prev() const;
    ///@}

public:

    Agent() = default;
    Agent(Agent && p);
    Agent(const Agent & p);
    Agent & operator=(const Agent & other_agent);
    ~Agent();

    /**
     * @name Add/Remove Virus/Tool
     *
     * Any of these is ultimately reflected at the end of the iteration. A
     * virus given to an agent takes it through the infection transition out
     * of its stage (see `NaturalHistory`); moving an agent otherwise is done
     * with `Model::fire()`.
     *
     * @param tool Tool to add
     * @param virus Virus to add
     */
    ///@{
    void add_tool(Model & model, const Tool & tool);
    void set_virus(Model & model, const Virus & virus);
    void rm_tool(Model & model, epiworld_fast_uint tool_idx);
    void rm_tool(Model & model, ToolPtr & tool);
    void rm_agent_by_virus(Model & model) = delete; ///< Agent removed by virus
    ///@}
    
    /**
     * @name Get the rates (multipliers) for the agent
     * 
     * @param v A pointer to a virus.
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_susceptibility_reduction(VirusPtr & v, Model & model);
    epiworld_double get_transmission_reduction(VirusPtr & v, Model & model);
    epiworld_double get_recovery_enhancer(VirusPtr & v, Model & model);
    epiworld_double get_death_reduction(VirusPtr & v, Model & model);
    ///@}

    int get_id() const; ///< Id of the individual

    VirusPtr & get_virus();
    const VirusPtr & get_virus() const;

    ToolPtr & get_tool(int i);
    ToolPtr & get_tool(std::string name);

    std::vector<ToolPtr> get_tools();
    const std::vector<ToolPtr> get_tools() const;
    size_t get_n_tools() const noexcept;

    /// @brief Applies the virus' mutation function (throws if the agent has no virus).
    void mutate_virus(Model & model);

    /// The last day the agent's stage or status changed.
    int get_state_last_changed() const;


    bool has_tool(epiworld_fast_uint t) const;
    bool has_tool(std::string_view name) const;
    bool has_tool(const Tool & t) const;
    bool has_virus(epiworld_fast_uint t) const;
    bool has_virus(std::string_view name) const;
    bool has_virus(const Virus & v) const;

    void print(Model & model, bool compressed = false) const;

    /**
     * @brief Access the j-th column of the agent
     * 
     * If an external array has been specified, then these two
     * functions can be used to access additional agent's features 
     * not included in the model.
     * 
     * @param j 
     * @param model Reference to the Model
     * @return double& 
     */
    ///@{
    double & operator()(size_t j, Model & model);
    double operator()(size_t j, const Model & model) const;
    ///@}

    bool operator==(const Agent & other) const;
    bool operator!=(const Agent & other) const {return !operator==(other);};

};



#endif
