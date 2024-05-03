#ifndef EPIWORLD_ENTITY_BONES_HPP
#define EPIWORLD_ENTITY_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Agent;

template<typename TSeq>
class AgentsSample;

template<typename TSeq>
inline void default_add_entity(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_entity(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
class Entity {
    friend class Agent<TSeq>;
    friend class AgentsSample<TSeq>;
    friend class Model<TSeq>;
    friend void default_add_entity<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_entity<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
private:

    Model<TSeq> * model;
    
    int id = -1;
    std::vector< size_t > agents;   ///< Vector of agents
    std::vector< size_t > agents_location; ///< Location where the entity is stored in the agent
    size_t n_agents = 0u;

    /**
     * @name Auxiliary variables for AgentsSample<TSeq> iterators
     * 
     * @details These variables+objects are used by the AgentsSample<TSeq>
     * class for building efficient iterators over agents. The idea is to
     * reduce the memory allocation, so only during the first call of
     * AgentsSample<TSeq>::AgentsSample(Entity<TSeq>) these vectors are allocated.
     */
    ///@{
    std::vector< Agent<TSeq> * > sampled_agents;
    size_t sampled_agents_n = 0u;
    std::vector< size_t > sampled_agents_left;
    size_t sampled_agents_left_n = 0u;
    // int date_last_add_or_remove = -99; ///< Last time the entity added or removed an agent
    ///@}

    int max_capacity = -1;
    std::string entity_name = "Unknown entity";

    std::vector< epiworld_double > location = {0.0}; ///< An arbitrary vector for location

    epiworld_fast_int state_init = -99;
    epiworld_fast_int state_post = -99;

    epiworld_fast_int queue_init = 0; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = 0; ///< Change of state when removed from agent.


public:

    epiworld_double prevalence = 0.0;
    bool prevalence_as_proportion = false;
    EntityToAgentFun<TSeq> dist_fun = nullptr;

    // Entity() = delete;
    // Entity(Entity<TSeq> & e) = delete;
    // Entity(const Entity<TSeq> & e);
    // Entity(Entity && e);
    Entity(std::string name) : entity_name(name) {};
    // Entity<TSeq> & operator=(const Entity<TSeq> & e);

    void add_agent(Agent<TSeq> & p, Model<TSeq> * model);
    void add_agent(Agent<TSeq> * p, Model<TSeq> * model);
    void rm_agent(size_t idx);
    size_t size() const noexcept;
    void set_location(std::vector< epiworld_double > loc);
    std::vector< epiworld_double > & get_location();

    typename std::vector< Agent<TSeq> * >::iterator begin();
    typename std::vector< Agent<TSeq> * >::iterator end();

    typename std::vector< Agent<TSeq> * >::const_iterator begin() const;
    typename std::vector< Agent<TSeq> * >::const_iterator end() const;

    Agent<TSeq> * operator[](size_t i);

    int get_id() const noexcept;
    const std::string & get_name() const noexcept;

    void set_state(epiworld_fast_int init, epiworld_fast_int post);
    void set_queue(epiworld_fast_int init, epiworld_fast_int post);
    void get_state(epiworld_fast_int * init, epiworld_fast_int * post);
    void get_queue(epiworld_fast_int * init, epiworld_fast_int * post);

    void reset();

    bool operator==(const Entity<TSeq> & other) const;
    bool operator!=(const Entity<TSeq> & other) const {return !operator==(other);};

    void distribute();

    std::vector< size_t > & get_agents();

    void print() const;

};


#endif
