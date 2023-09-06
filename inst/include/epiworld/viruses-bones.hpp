#ifndef EPIWORLD_VIRUSES_BONES_HPP
#define EPIWORLD_VIRUSES_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;

/**
 * @brief Set of viruses (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Viruses {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< VirusPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_viruses;

public:

    Viruses() = delete;
    Viruses(Agent<TSeq> & p) : dat(&p.viruses), n_viruses(&p.n_viruses) {};

    typename std::vector< VirusPtr<TSeq> >::iterator begin();
    typename std::vector< VirusPtr<TSeq> >::iterator end();

    VirusPtr<TSeq> & operator()(size_t i);
    VirusPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    void print() const noexcept;

};

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::iterator Viruses<TSeq>::begin()
{

    if (*n_viruses == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::iterator Viruses<TSeq>::end()
{
    
    #ifdef EPI_DEBUG
    if (dat->size() < *n_viruses)
        throw EPI_DEBUG_ERROR(std::logic_error, "Viruses:: The end of the virus is out of range");
    #endif 

    return begin() + *n_viruses;
}

template<typename TSeq>
inline VirusPtr<TSeq> & Viruses<TSeq>::operator()(size_t i)
{

    if (i >= *n_viruses)
        throw std::range_error("Virus index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline VirusPtr<TSeq> & Viruses<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Viruses<TSeq>::size() const noexcept 
{
    return *n_viruses;
}

template<typename TSeq>
inline void Viruses<TSeq>::print() const noexcept
{

    if (*n_viruses == 0u)
    {
        printf_epiworld("List of viruses (none)\n");
        return;
    }

    printf_epiworld("List of viruses (%i): ", *n_viruses);

    // Printing the name of each virus separated by a comma
    for (size_t i = 0u; i < *n_viruses; ++i)
    {
        if (i == *n_viruses - 1u)
        {
            printf_epiworld("%s", dat->operator[](i)->get_name().c_str());
        } else 
        {
            printf_epiworld("%s, ", dat->operator[](i)->get_name().c_str());
        }
    }
    
    printf_epiworld("\n");

}

/**
 * @brief Set of Viruses (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Viruses_const {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< VirusPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_viruses;

public:

    Viruses_const() = delete;
    Viruses_const(const Agent<TSeq> & p) : dat(&p.viruses), n_viruses(&p.n_viruses) {};

    typename std::vector< VirusPtr<TSeq> >::const_iterator begin() const;
    typename std::vector< VirusPtr<TSeq> >::const_iterator end() const;

    const VirusPtr<TSeq> & operator()(size_t i);
    const VirusPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    void print() const noexcept;

};

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::const_iterator Viruses_const<TSeq>::begin() const {

    if (*n_viruses == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< VirusPtr<TSeq> >::const_iterator Viruses_const<TSeq>::end() const {

    #ifdef EPI_DEBUG
    if (dat->size() < *n_viruses)
        throw EPI_DEBUG_ERROR(std::logic_error, "Viruses_const:: The end of the virus is out of range");
    #endif 
    return begin() + *n_viruses;
}

template<typename TSeq>
inline const VirusPtr<TSeq> & Viruses_const<TSeq>::operator()(size_t i)
{

    if (i >= *n_viruses)
        throw std::range_error("Virus index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline const VirusPtr<TSeq> & Viruses_const<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Viruses_const<TSeq>::size() const noexcept 
{
    return *n_viruses;
}

template<typename TSeq>
inline void Viruses_const<TSeq>::print() const noexcept
{

    if (*n_viruses == 0u)
    {
        printf_epiworld("List of viruses (none)\n");
        return;
    }

    printf_epiworld("List of viruses (%i): ", *n_viruses);

    // Printing the name of each virus separated by a comma
    for (size_t i = 0u; i < *n_viruses; ++i)
    {
        if (i == *n_viruses - 1u)
        {
            printf_epiworld("%s", dat->operator[](i)->get_name().c_str());
        } else
        {
            printf_epiworld("%s, ", dat->operator[](i)->get_name().c_str());
        }
    }
    
    printf_epiworld("\n");

}


#endif