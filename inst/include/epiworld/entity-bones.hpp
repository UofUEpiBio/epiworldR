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
    
    int id = -1;
    std::vector< size_t > agents;   ///< Vector of agents
    std::vector< size_t > agents_location; ///< Location where the entity is stored in the agent
    size_t n_agents = 0u;

    int max_capacity = -1;
    std::string entity_name = "Unnamed entity";

    std::vector< epiworld_double > location = {0.0}; ///< An arbitrary vector for location

    epiworld_fast_int state_init = -99;
    epiworld_fast_int state_post = -99;

    epiworld_fast_int queue_init = 0; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = 0; ///< Change of state when removed from agent.

    EntityToAgentFun<TSeq> dist_fun = nullptr;

public:


    /**
     * @brief Constructs an Entity object.
     *
     * This constructor initializes an Entity object with the specified parameters.
     *
     * @param name The name of the entity.
     * @param fun A function pointer to a function that maps the entity to an agent.
     */
    Entity(
        std::string name,
        EntityToAgentFun<TSeq> fun = nullptr
        ) :
            entity_name(name),
            dist_fun(fun)
        {};
    
    void add_agent(Agent<TSeq> & p, Model<TSeq> * model);
    void add_agent(Agent<TSeq> * p, Model<TSeq> * model);
    void rm_agent(size_t idx, Model<TSeq> * model);
    size_t size() const noexcept;
    void set_location(std::vector< epiworld_double > loc);
    std::vector< epiworld_double > & get_location();

    typename std::vector< Agent<TSeq> * >::iterator begin();
    typename std::vector< Agent<TSeq> * >::iterator end();

    typename std::vector< Agent<TSeq> * >::const_iterator begin() const;
    typename std::vector< Agent<TSeq> * >::const_iterator end() const;

    size_t operator[](size_t i);

    int get_id() const noexcept;
    const std::string & get_name() const noexcept;

    void set_state(epiworld_fast_int init, epiworld_fast_int post);
    void set_queue(epiworld_fast_int init, epiworld_fast_int post);
    void get_state(epiworld_fast_int * init, epiworld_fast_int * post);
    void get_queue(epiworld_fast_int * init, epiworld_fast_int * post);

    void reset();

    bool operator==(const Entity<TSeq> & other) const;
    bool operator!=(const Entity<TSeq> & other) const {return !operator==(other);};

    /** 
     * @name Entity distribution
     * 
     * @details These functions are used for distributing agents among entities.
     * The idea is to have a flexible way of distributing agents among entities.
     
     */
    void distribute(Model<TSeq> * model);

    std::vector< size_t > & get_agents();

    void print() const;
    void set_distribution(EntityToAgentFun<TSeq> fun);

};


#endif
