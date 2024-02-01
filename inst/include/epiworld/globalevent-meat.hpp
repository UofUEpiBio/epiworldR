#ifndef EPIWORLD_GLOBALEVENT_MEAT_HPP
#define EPIWORLD_GLOBALEVENT_MEAT_HPP

template<typename TSeq>
inline GlobalEvent<TSeq>::GlobalEvent(
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
inline void GlobalEvent<TSeq>::operator()(Model<TSeq> * m, int day)
{   
    
    if (this->fun == nullptr)
        return;

    // events apply if day is negative or if day is equal to the day of the action
    if (this->day < 0 || this->day == day)
        this->fun(m);
    
    return;

}

template<typename TSeq>
inline void GlobalEvent<TSeq>::set_name(std::string name)
{
    this->name = name;
}

template<typename TSeq>
inline std::string GlobalEvent<TSeq>::get_name() const
{
    return this->name;
}

template<typename TSeq>
inline void GlobalEvent<TSeq>::set_day(int day)
{
    this->day = day;
}

template<typename TSeq>
inline int GlobalEvent<TSeq>::get_day() const
{
    return this->day;
}

template<typename TSeq>
inline void GlobalEvent<TSeq>::print() const
{
    printf_epiworld(
        "Global action: %s\n"
        "  - Day: %i\n",
        this->name.c_str(),
        this->day
        );
}

template<typename TSeq>
inline bool GlobalEvent<TSeq>::operator==(const GlobalEvent<TSeq> & other) const
{
    return (this->name == other.name) && (this->day == other.day);
}

template<typename TSeq>
inline bool GlobalEvent<TSeq>::operator!=(const GlobalEvent<TSeq> & other) const
{
    return !(*this == other);
}

#endif
