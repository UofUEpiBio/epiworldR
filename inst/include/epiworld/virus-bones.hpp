#ifndef EPIWORLD_VIRUS_HPP
#define EPIWORLD_VIRUS_HPP

template<typename TSeq>
class Agent;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Model;

/**
 * @brief Virus
 * 
 * @tparam TSeq 
 * @details
 * Raw transmisibility of a virus should be a function of its genetic
 * sequence. Nonetheless, transmisibility can be reduced as a result of
 * having one or more tools to fight the virus. Because of this, transmisibility
 * should be a function of the agent.
 */
template<typename TSeq>
class Virus {
    friend class Agent<TSeq>;
    friend class Model<TSeq>;
    friend class DataBase<TSeq>;
    friend void default_add_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Action<TSeq> & a, Model<TSeq> * m);
private:
    
    Agent<TSeq> * agent       = nullptr;
    int       pos_in_agent    = -99; ///< Location in the agent
    int agent_exposure_number = -99;

    std::shared_ptr<TSeq> baseline_sequence = nullptr;
    std::shared_ptr<std::string> virus_name = nullptr;
    int date = -99;
    int id   = -99;
    bool active = true;
    MutFun<TSeq>          mutation_fun                 = nullptr;
    PostRecoveryFun<TSeq> post_recovery_fun            = nullptr;
    VirusFun<TSeq>        probability_of_infecting_fun = nullptr;
    VirusFun<TSeq>        probability_of_recovery_fun  = nullptr;
    VirusFun<TSeq>        probability_of_death_fun     = nullptr;
    VirusFun<TSeq>        incubation_fun               = nullptr;

    // Setup parameters
    std::vector< epiworld_double * > params = {};
    std::vector< epiworld_double > data = {};

    epiworld_fast_int state_init    = -99; ///< Change of state when added to agent.
    epiworld_fast_int state_post    = -99; ///< Change of state when removed from agent.
    epiworld_fast_int state_removed = -99; ///< Change of state when agent is removed

    epiworld_fast_int queue_init    = Queue<TSeq>::Everyone; ///< Change of state when added to agent.
    epiworld_fast_int queue_post    = -Queue<TSeq>::Everyone; ///< Change of state when removed from agent.
    epiworld_fast_int queue_removed = -99; ///< Change of state when agent is removed

public:
    Virus(std::string name = "unknown virus");

    void mutate(Model<TSeq> * model);
    void set_mutation(MutFun<TSeq> fun);
    
    std::shared_ptr<TSeq> get_sequence();
    void set_sequence(TSeq sequence);
    
    Agent<TSeq> * get_agent();
    void set_agent(Agent<TSeq> * p, epiworld_fast_uint idx);
    
    void set_date(int d);
    int get_date() const;

    void set_id(int idx);
    int get_id() const;

    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_prob_infecting(Model<TSeq> * model);
    epiworld_double get_prob_recovery(Model<TSeq> * model);
    epiworld_double get_prob_death(Model<TSeq> * model);
    epiworld_double get_incubation(Model<TSeq> * model);
    
    void post_recovery(Model<TSeq> * model);
    void set_post_recovery(PostRecoveryFun<TSeq> fun);
    void set_post_immunity(epiworld_double prob);
    void set_post_immunity(epiworld_double * prob);

    void set_prob_infecting_fun(VirusFun<TSeq> fun);
    void set_prob_recovery_fun(VirusFun<TSeq> fun);
    void set_prob_death_fun(VirusFun<TSeq> fun);
    void set_incubation_fun(VirusFun<TSeq> fun);
    
    void set_prob_infecting(const epiworld_double * prob);
    void set_prob_recovery(const epiworld_double * prob);
    void set_prob_death(const epiworld_double * prob);
    void set_incubation(const epiworld_double * prob);
    
    void set_prob_infecting(epiworld_double prob);
    void set_prob_recovery(epiworld_double prob);
    void set_prob_death(epiworld_double prob);
    void set_incubation(epiworld_double prob);
    ///@}


    void set_name(std::string name);
    std::string get_name() const;

    std::vector< epiworld_double > & get_data();

    /**
     * @name Get and set the state and queue
     * 
     * After applied, viruses can change the state and affect
     * the queue of agents. These function sets the default values,
     * which are retrieved when adding or removing a virus does not
     * specify a change in state or in queue.
     * 
     * @param init After the virus/tool is added to the agent.
     * @param end After the virus/tool is removed.
     * @param removed After the agent (Agent) is removed.
     */
    ///@{
    void set_state(
        epiworld_fast_int init,
        epiworld_fast_int end,
        epiworld_fast_int removed = -99
        );
        
    void set_queue(
        epiworld_fast_int init,
        epiworld_fast_int end,
        epiworld_fast_int removed = -99
        );

    void get_state(
        epiworld_fast_int * init,
        epiworld_fast_int * end,
        epiworld_fast_int * removed = nullptr
        );

    void get_queue(
        epiworld_fast_int * init,
        epiworld_fast_int * end,
        epiworld_fast_int * removed = nullptr
        );
    ///@}

    bool operator==(const Virus<TSeq> & other) const;
    bool operator!=(const Virus<TSeq> & other) const {return !operator==(other);};

    void print() const;

};

#endif