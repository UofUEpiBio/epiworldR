#ifndef EPIWORLD_TOOLS_BONES_HPP
#define EPIWORLD_TOOLS_BONES_HPP

template<typename TSeq>
class Tool;

template<typename TSeq>
class Agent;

// #define ToolPtr<TSeq> std::shared_ptr< Tool<TSeq> >

/**
 * @brief Set of tools (useful for building iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Tools {
    friend class Tool<TSeq>;
    friend class Agent<TSeq>;
private:
    std::vector< ToolPtr<TSeq> > * dat;
    const unsigned int * n_tools;

public:

    Tools() = delete;
    Tools(Agent<TSeq> & p) : dat(&p.tools), n_tools(&p.n_tools) {};

    typename std::vector< ToolPtr<TSeq> >::iterator begin();
    typename std::vector< ToolPtr<TSeq> >::iterator end();

    ToolPtr<TSeq> & operator()(size_t i);
    ToolPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    void print() const noexcept;

};

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::iterator Tools<TSeq>::begin()
{

    if (*n_tools == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::iterator Tools<TSeq>::end()
{
     
    return begin() + *n_tools;
}

template<typename TSeq>
inline ToolPtr<TSeq> & Tools<TSeq>::operator()(size_t i)
{

    if (i >= *n_tools)
        throw std::range_error("Tool index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline ToolPtr<TSeq> & Tools<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Tools<TSeq>::size() const noexcept 
{
    return *n_tools;
}

template<typename TSeq>
inline void Tools<TSeq>::print() const noexcept 
{
    if (*n_tools == 0u)
    {
        printf_epiworld("List of tools (none)\n");
        return;
    }

    printf_epiworld("List of tools (%i): ", static_cast<int>(*n_tools));

    // Printing the name of each virus separated by a comma
    for (size_t i = 0u; i < *n_tools; ++i)
    {
        if (i == *n_tools - 1u)
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
 * @brief Set of Tools (const) (useful for iterators)
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Tools_const {
    friend class Tool<TSeq>;
    friend class Agent<TSeq>;
private:
    const std::vector< ToolPtr<TSeq> > * dat;
    const epiworld_fast_uint * n_tools;

public:

    Tools_const() = delete;
    Tools_const(const Agent<TSeq> & p) : dat(&p.tools), n_tools(&p.n_tools) {};

    typename std::vector< ToolPtr<TSeq> >::const_iterator begin() const;
    typename std::vector< ToolPtr<TSeq> >::const_iterator end() const;

    const ToolPtr<TSeq> & operator()(size_t i);
    const ToolPtr<TSeq> & operator[](size_t i);

    size_t size() const noexcept;

    void print() const noexcept;

};

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::const_iterator Tools_const<TSeq>::begin() const {

    if (*n_tools == 0u)
        return dat->end();
    
    return dat->begin();
}

template<typename TSeq>
inline typename std::vector< ToolPtr<TSeq> >::const_iterator Tools_const<TSeq>::end() const {
     
    return begin() + *n_tools;
}

template<typename TSeq>
inline const ToolPtr<TSeq> & Tools_const<TSeq>::operator()(size_t i)
{

    if (i >= *n_tools)
        throw std::range_error("Tool index out of range.");

    return dat->operator[](i);

}

template<typename TSeq>
inline const ToolPtr<TSeq> & Tools_const<TSeq>::operator[](size_t i)
{

    return dat->operator[](i);

}

template<typename TSeq>
inline size_t Tools_const<TSeq>::size() const noexcept 
{
    return *n_tools;
}

template<typename TSeq>
inline void Tools_const<TSeq>::print() const noexcept 
{
    if (*n_tools == 0u)
    {
        printf_epiworld("List of tools (none)\n");
        return;
    }

    printf_epiworld("List of tools (%i): ", *n_tools);

    // Printing the name of each virus separated by a comma
    for (size_t i = 0u; i < *n_tools; ++i)
    {
        if (i == *n_tools - 1u)
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