#ifndef EPIWORLD_GLOBALACTIONS_BONES_HPP
#define EPIWORLD_GLOBALACTIONS_BONES_HPP

// template<typename TSeq = EPI_DEFAULT_TSEQ>
// using GlobalFun = std::function<void(Model<TSeq>*)>;

/**
 * @brief Template for a Global Action
 * @details Global actions are functions that Model<TSeq> executes
 * at the end of a day.
 * 
 */
template<typename TSeq>
class GlobalAction
{
private:
    GlobalFun<TSeq> fun = nullptr;
    std::string name = "A global action";
    int day = -99;
public:

    GlobalAction() {};

    /**
     * @brief Construct a new Global Action object
     * 
     * @param fun A function that takes a Model<TSeq> * as argument and returns void.
     * @param name A descriptive name for the action.
     * @param day The day when the action will be executed. If negative, it will be executed every day.
     */
    GlobalAction(GlobalFun<TSeq> fun, std::string name, int day = -99);
    
    ~GlobalAction() {};

    void operator()(Model<TSeq> * m, int day);

    void set_name(std::string name);
    std::string get_name() const;

    void set_day(int day);
    int get_day() const;
    
    void print() const;

    // Comparison operators
    bool operator==(const GlobalAction<TSeq> & other) const;
    bool operator!=(const GlobalAction<TSeq> & other) const;

};



#endif
