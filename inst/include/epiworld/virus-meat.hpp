#ifndef EPIWORLD_VIRUS_MEAT_HPP
#define EPIWORLD_VIRUS_MEAT_HPP

#include "config.hpp"

inline Virus::Virus()
{

    #ifdef EPI_DEBUG_VIRUS
    counter_construct++;
    #endif

}

inline Virus::Virus(
    std::string name
    ) {

    #ifdef EPI_DEBUG_VIRUS
    counter_construct++;
    #endif
    
    set_name(name);
}

#ifdef EPI_DEBUG_VIRUS
inline std::atomic<int> Virus::counter_construct = 0;

inline std::atomic<int> Virus::counter_copy_construct = 0;

inline std::atomic<int> Virus::counter_move_construct = 0;

inline std::atomic<int> Virus::counter_copy_assign = 0;

inline std::atomic<int> Virus::counter_move_assign = 0;

inline std::atomic<int> Virus::counter_destruct = 0;

inline Virus::~Virus()
{
    counter_destruct++;
}

// Copy constructor
inline Virus::Virus(const Virus& other)
    : agent(other.agent),
      virus_name(other.virus_name),
      date(other.date),
      id(other.id),
      lineage_id(other.lineage_id),
      lineage_bit(other.lineage_bit),
      param_redirects(other.param_redirects),
      nh_def(other.nh_def)
{
    counter_copy_construct++;
}

// Move constructor
inline Virus::Virus(Virus&& other) noexcept
    : agent(other.agent),
      virus_name(std::move(other.virus_name)),
      date(other.date),
      id(other.id),
      lineage_id(other.lineage_id),
      lineage_bit(other.lineage_bit),
      param_redirects(std::move(other.param_redirects)),
      nh_def(std::move(other.nh_def))
{
    counter_move_construct++;
    // other.agent = nullptr;
}

// Copy assignment
inline Virus& Virus::operator=(const Virus& other)
{
    if (this != &other) {
        agent = other.agent;
        virus_name = other.virus_name;
        date = other.date;
        id = other.id;
        lineage_id = other.lineage_id;
        lineage_bit = other.lineage_bit;
        param_redirects = other.param_redirects;
        nh_def = other.nh_def;
        counter_copy_assign++;
    }
    return *this;
}

// Move assignment
inline Virus& Virus::operator=(Virus&& other) noexcept
{
    if (this != &other) {
        agent = other.agent;
        virus_name = std::move(other.virus_name);
        date = other.date;
        id = other.id;
        lineage_id = other.lineage_id;
        lineage_bit = other.lineage_bit;
        param_redirects = std::move(other.param_redirects);
        nh_def = std::move(other.nh_def);
        other.agent = nullptr;
        counter_move_assign++;
    }
    return *this;
}
#endif

inline void Virus::mutate(
    Model * model
) {

    if (mutation)
        if (mutation(agent, *this, model))
        {
            check_params(*model);
            model->get_db().record_virus(*this);
        }

    return;
    
}

inline void Virus::set_mutation(
    MutFun fun
) {
    mutation = MutFun(fun);
}


inline Agent * Virus::get_agent()
{

    return agent;

}

inline void Virus::set_agent(Agent * p)
{
    agent = p;
}

inline void Virus::set_id(int idx)
{

    id = idx;
    return;

}

inline int Virus::get_id() const
{
    
    return id;

}

inline int Virus::get_lineage_id() const
{
    return lineage_id;
}

inline void Virus::set_date(int d) 
{

    date = d;
    return;

}

inline int Virus::get_date() const
{
    
    return date;
    
}

inline void Virus::set_name(std::string name)
{

    virus_name = name;

}

inline std::string Virus::get_name() const
{

    return virus_name;

}

inline bool Virus::operator==(const Virus & other) const
{
    
    EPI_DEBUG_FAIL_AT_TRUE(
        virus_name != other.virus_name,
        "Virus:: virus_name don't match"
    )
    






    return true;

}

inline void Virus::print() const
{

    printf_epiworld("Virus         : %s\n", virus_name.c_str());
    printf_epiworld("Id            : %s\n", (id < 0)? std::string("(empty)").c_str() : std::to_string(id).c_str());

}

inline std::unique_ptr<Virus> Virus::clone_ptr() const
{
    auto cloned = std::make_unique<Virus>(*this);
    return cloned;
}

inline void Virus::use_param(std::string param, std::string replacement)
{
    for (auto & r : param_redirects)
        if (r.first.name() == param)
        {
            r.second = ParamRef(std::move(replacement));
            return;
        }

    param_redirects.emplace_back(
        ParamRef(std::move(param)), ParamRef(std::move(replacement))
    );
}

inline ParamId Virus::param_redirect(
    ParamId id,
    const Model & model
) const
{
    for (const auto & r : param_redirects)
        if (r.first.id(model).idx == id.idx)
            return r.second.id(model);

    return id;
}

inline NaturalHistory & Virus::nh_edit()
{
    // A virus in a model (its copy there, or an agent's) keeps the model's
    // stages and transitions; a strain changes parameters instead
    if (id >= 0)
        throw std::logic_error(
            "The virus \"" + virus_name + "\" is in a model: its natural "
            "history cannot change (a strain reads other parameters with "
            "use_param())."
        );

    // Copies share the natural history until one of them is edited
    if (!nh_def)
        nh_def = std::make_shared< NaturalHistory >(virus_name);
    else if (nh_def.use_count() > 1)
        nh_def = std::make_shared< NaturalHistory >(*nh_def);

    return *nh_def;
}

inline int Virus::add_stage(std::string label, unsigned int flags)
{
    return nh_edit().add_stage(std::move(label), flags);
}

inline void Virus::set_contact(int s, epiworld_double weight)
{
    nh_edit().set_contact(s, weight);
}

inline TransitionRef Virus::infection(
    Nodes from, int to, Value per_contact
)
{
    return nh_edit().infection(std::move(from), to, std::move(per_contact));
}

inline TransitionRef Virus::rate(
    Nodes from, int to, Value daily_probability
)
{
    return nh_edit().rate(std::move(from), to, std::move(daily_probability));
}

inline TransitionRef Virus::after(int from, int to, Duration d)
{
    return nh_edit().after(from, to, std::move(d));
}

inline TransitionRef Virus::external(
    Nodes from, int to, std::string label
)
{
    return nh_edit().external(std::move(from), to, std::move(label));
}

inline void Virus::check_params(const Model & model) const
{
    for (const auto & r : param_redirects)
        if (!model.has_param(r.second.name()))
            throw std::logic_error(
                "The virus \"" + virus_name + "\" reads the parameter \"" +
                r.second.name() + "\" in place of \"" + r.first.name() +
                "\", but the model has no such parameter (add it with "
                "add_param())."
            );
}

#endif
