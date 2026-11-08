#ifndef EPIWORLD_GLOBALEVENT_MEAT_HPP
#define EPIWORLD_GLOBALEVENT_MEAT_HPP

#include <memory>
#include "globalevent-bones.hpp"

inline GlobalEvent::GlobalEvent(
    GlobalFun fun,
    std::string name,
    int day
    )
{
    this->fun = fun;
    this->name = name;
    this->day = day;
}

inline void GlobalEvent::operator()(Model * m, int day)
{   
    
    if (this->fun == nullptr)
        return;

    // events apply if day is negative or if day is equal to the day of the action
    if (this->day < 0 || this->day == day)
        this->fun(m);
    
    return;

}

inline void GlobalEvent::reset(Model *)
{
    return;
}

inline void GlobalEvent::set_name(std::string name)
{
    this->name = name;
}

inline std::string GlobalEvent::get_name() const
{
    return this->name;
}

inline void GlobalEvent::set_day(int day)
{
    this->day = day;
}

inline int GlobalEvent::get_day() const
{
    return this->day;
}

inline void GlobalEvent::print() const
{
    printf_epiworld(
        "Global action: %s\n"
        "  - Day: %i\n",
        this->name.c_str(),
        this->day
        );
}

inline bool GlobalEvent::operator==(const GlobalEvent & other) const
{
    return (this->name == other.name) && (this->day == other.day);
}

inline bool GlobalEvent::operator!=(const GlobalEvent & other) const
{
    return !(*this == other);
}

inline std::unique_ptr<GlobalEvent> GlobalEvent::clone_ptr() const
{
    return std::make_unique<GlobalEvent>(*this);
}

#endif
