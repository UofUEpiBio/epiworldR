#ifndef EPIWORLD_GLOBALACTIONS_MEAT_HPP
#define EPIWORLD_GLOBALACTIONS_MEAT_HPP

template<typename TSeq>
inline GlobalAction<TSeq>::GlobalAction(
    GlobalFun<TSeq> fun,
    std::string name,
    int day
    )
{
    this->fun = fun;
    this->name = name;
    this->day = day;
}

template<typename TSeq>
inline void GlobalAction<TSeq>::operator()(Model<TSeq> * m, int day)
{   
    
    if (this->fun == nullptr)
        return;

    // Actions apply if day is negative or if day is equal to the day of the action
    if (this->day < 0 || this->day == day)
        this->fun(m);
    
    return;

}

template<typename TSeq>
inline void GlobalAction<TSeq>::set_name(std::string name)
{
    this->name = name;
}

template<typename TSeq>
inline std::string GlobalAction<TSeq>::get_name() const
{
    return this->name;
}

template<typename TSeq>
inline void GlobalAction<TSeq>::set_day(int day)
{
    this->day = day;
}

template<typename TSeq>
inline int GlobalAction<TSeq>::get_day() const
{
    return this->day;
}

template<typename TSeq>
inline void GlobalAction<TSeq>::print() const
{
    printf_epiworld(
        "Global action: %s\n"
        "  - Day: %i\n",
        this->name.c_str(),
        this->day
        );
}

template<typename TSeq>
inline bool GlobalAction<TSeq>::operator==(const GlobalAction<TSeq> & other) const
{
    return (this->name == other.name) && (this->day == other.day);
}

template<typename TSeq>
inline bool GlobalAction<TSeq>::operator!=(const GlobalAction<TSeq> & other) const
{
    return !(*this == other);
}

#endif
