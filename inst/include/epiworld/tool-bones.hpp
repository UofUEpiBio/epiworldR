
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
    friend void default_add_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
private:

    Agent<TSeq> * agent = nullptr;
    int pos_in_agent        = -99; ///< Location in the agent

    int date = -99;
    int id   = -99;
    std::shared_ptr<std::string> tool_name     = nullptr;
    std::shared_ptr<TSeq> sequence             = nullptr;
    ToolFun<TSeq> susceptibility_reduction_fun = nullptr;
    ToolFun<TSeq> transmission_reduction_fun   = nullptr;
    ToolFun<TSeq> recovery_enhancer_fun        = nullptr;
    ToolFun<TSeq> death_reduction_fun          = nullptr;

    // Setup parameters
    std::vector< epiworld_double * > params;  

    epiworld_fast_int state_init = -99;
    epiworld_fast_int state_post = -99;

    epiworld_fast_int queue_init = Queue<TSeq>::NoOne; ///< Change of state when added to agent.
    epiworld_fast_int queue_post = Queue<TSeq>::NoOne; ///< Change of state when removed from agent.

    void set_agent(Agent<TSeq> * p, size_t idx);

    ToolToAgentFun<TSeq> dist_fun = nullptr;

public:
    Tool(std::string name = "unknown tool");
    Tool(
        std::string name,
        epiworld_double prevalence,
        bool as_proportion
    );

    void set_sequence(TSeq d);
    void set_sequence(std::shared_ptr<TSeq> d);
    std::shared_ptr<TSeq> get_sequence();

    /**
     * @name Get and set the tool functions
     * 
     * @param v The virus over which to operate
     * @param fun the function to be used
     * 
     * @return epiworld_double 
     */
    ///@{
    epiworld_double get_susceptibility_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_transmission_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_recovery_enhancer(VirusPtr<TSeq> v, Model<TSeq> * model);
    epiworld_double get_death_reduction(VirusPtr<TSeq> v, Model<TSeq> * model);
    
    void set_susceptibility_reduction_fun(ToolFun<TSeq> fun);
    void set_transmission_reduction_fun(ToolFun<TSeq> fun);
    void set_recovery_enhancer_fun(ToolFun<TSeq> fun);
    void set_death_reduction_fun(ToolFun<TSeq> fun);

    void set_susceptibility_reduction(epiworld_double * prob);
    void set_transmission_reduction(epiworld_double * prob);
    void set_recovery_enhancer(epiworld_double * prob);
    void set_death_reduction(epiworld_double * prob);

    void set_susceptibility_reduction(epiworld_double prob);
    void set_transmission_reduction(epiworld_double prob);
    void set_recovery_enhancer(epiworld_double prob);
    void set_death_reduction(epiworld_double prob);
    ///@}

    void set_name(std::string name);
    std::string get_name() const;

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

};

#endif