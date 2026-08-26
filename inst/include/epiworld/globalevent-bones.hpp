#ifndef EPIWORLD_GLOBALEVENT_BONES_HPP
#define EPIWORLD_GLOBALEVENT_BONES_HPP

#include <string>
#include "config.hpp"

/**
 * @brief Template for a Global Event
 * @details Global events are functions that Model<TSeq> executes
 * at the end of a day.
 * 
 */
template<typename TSeq>
class GlobalEvent
{
private:
    GlobalFun<TSeq> fun = nullptr;
    std::string name = "A global action";
    int day = -99;
public:

    GlobalEvent() {};

    /**
     * @brief Construct a new Global Event object
     * 
     * @param fun A function that takes a Model<TSeq> * as argument and returns void.
     * @param name A descriptive name for the action.
     * @param day The day when the action will be executed. If negative, it will be executed every day.
     */
    GlobalEvent(GlobalFun<TSeq> fun, std::string name, int day = -99);
    
    virtual ~GlobalEvent() = default;

    virtual void operator()(Model<TSeq> * m, int day);

    /**
     * @brief Prepare the event for a run.
     *
     * @details Called by `Model::reset()` at the end of every reset, i.e. once
     * per run and once per replicate of `Model::run_multiple()`, with the run's
     * RNG already seeded and the agents, tools, and viruses distributed. Events
     * that need to set themselves up on the model they are running in -- the
     * moment before day 1 -- do it here; the default does nothing.
     *
     * @param m The model the event belongs to.
     */
    virtual void reset(Model<TSeq> * m);

    void set_name(std::string name);
    std::string get_name() const;

    void set_day(int day);
    int get_day() const;
    
    void print() const;

    // Comparison operators
    bool operator==(const GlobalEvent<TSeq> & other) const;
    bool operator!=(const GlobalEvent<TSeq> & other) const;

    virtual std::unique_ptr<GlobalEvent<TSeq>> clone_ptr() const;

};



#endif
