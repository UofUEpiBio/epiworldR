#ifndef EPIWORLD_GLOBALEVENT_BONES_HPP
#define EPIWORLD_GLOBALEVENT_BONES_HPP

// template<typename TSeq = EPI_DEFAULT_TSEQ>
// using GlobalFun = std::function<void(Model<TSeq>*)>;

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
    
    ~GlobalEvent() {};

    void operator()(Model<TSeq> * m, int day);

    void set_name(std::string name);
    std::string get_name() const;

    void set_day(int day);
    int get_day() const;
    
    void print() const;

    // Comparison operators
    bool operator==(const GlobalEvent<TSeq> & other) const;
    bool operator!=(const GlobalEvent<TSeq> & other) const;

};



#endif
