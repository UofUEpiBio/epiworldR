#ifndef EPIWORLD_ENTITIES_BONES_HPP
#define EPIWORLD_ENTITIES_BONES_HPP

template<typename TSeq>
class Virus;

template<typename TSeq>
class Agent;


/**
 * @brief Set of Entities (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Entities {
    friend class Entity<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< Entity<TSeq> * >  dat;
    const size_t n_entities;

public:

    Entities() = delete;
    Entities(Agent<TSeq> & p);

    typename std::vector< Entity<TSeq> * >::iterator begin();
    typename std::vector< Entity<TSeq> * >::iterator end();

    Entity<TSeq> & operator()(size_t i);
    Entity<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    bool operator==(const Entities<TSeq> & other) const;

};

template<typename TSeq>
inline Entities<TSeq>::Entities(Agent<TSeq> & p) :
    n_entities(p.get_n_entities())
{

    dat.reserve(n_entities);
    for (size_t i = 0u; i < n_entities; ++i)
        dat.push_back(&p.get_entity(i));

}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::iterator Entities<TSeq>::begin()
{

    if (n_entities == 0u)
        return dat.end();
    
    return dat.begin();
}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::iterator Entities<TSeq>::end()
{
     
    return begin() + n_entities;
}

template<typename TSeq>
inline Entity<TSeq> & Entities<TSeq>::operator()(size_t i)
{

    if (i >= n_entities)
        throw std::range_error("Entity index out of range.");

    return *dat[i];

}

template<typename TSeq>
inline Entity<TSeq> & Entities<TSeq>::operator[](size_t i)
{

    return *dat[i];

}

template<typename TSeq>
inline size_t Entities<TSeq>::size() const noexcept 
{
    return n_entities;
}

template<typename TSeq>
inline bool Entities<TSeq>::operator==(const Entities<TSeq> & other) const
{

    if (n_entities != other.n_entities)
        return false;

    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (dat[i] != other.dat[i])
            return false;
    }

    return true;
}

/**
 * @brief Set of Entities (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Entities_const {
    friend class Virus<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< Entity<TSeq>* > dat;
    const size_t n_entities;

public:

    Entities_const() = delete;
    Entities_const(const Agent<TSeq> & p);

    typename std::vector< Entity<TSeq>* >::const_iterator begin();
    typename std::vector< Entity<TSeq>* >::const_iterator end();

    const Entity<TSeq> & operator()(size_t i);
    const Entity<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    bool operator==(const Entities_const<TSeq> & other) const;

};

template<typename TSeq>
inline Entities_const<TSeq>::Entities_const(const Agent<TSeq> & p) :
    n_entities(p.get_n_entities())
{

    dat.reserve(n_entities);
    for (size_t i = 0u; i < n_entities; ++i)
        dat.push_back(&p.get_entity(i));

}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::const_iterator Entities_const<TSeq>::begin() {

    if (n_entities == 0u)
        return dat.end();
    
    return dat.begin();
}

template<typename TSeq>
inline typename std::vector< Entity<TSeq>* >::const_iterator Entities_const<TSeq>::end() {
     
    return begin() + n_entities;
}

template<typename TSeq>
inline const Entity<TSeq> & Entities_const<TSeq>::operator()(size_t i)
{

    if (i >= n_entities)
        throw std::range_error("Entity index out of range.");

    return *dat[i];

}

template<typename TSeq>
inline const Entity<TSeq> & Entities_const<TSeq>::operator[](size_t i)
{

    return *dat[i];

}

template<typename TSeq>
inline size_t Entities_const<TSeq>::size() const noexcept 
{
    return n_entities;
}

template<typename TSeq>
inline bool Entities_const<TSeq>::operator==(const Entities_const<TSeq> & other) const
{
    
    if (n_entities != other.n_entities)
        return false;

    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (dat[i] != other.dat[i])
            return false;
    }

    return true;
}


#endif