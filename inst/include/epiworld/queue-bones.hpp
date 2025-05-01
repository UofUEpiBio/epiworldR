#ifndef EPIWORLD_QUEUE_BONES_HPP
#define EPIWORLD_QUEUE_BONES_HPP

/**
 * @brief Controls which agents are verified at each step
 * 
 * @details The idea is that only agents who are either in
 * an infected state or have an infected neighbor should be
 * checked. Otherwise it makes no sense (no chance to recover
 * or capture the disease).
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class Queue
{
    friend class Model<TSeq>;

private:

    /**
     * @brief Count of ego's neighbors in queue (including ego)
     */
    std::vector< epiworld_fast_int > active;
    Model<TSeq> * model = nullptr;
    int n_in_queue = 0;

    // Auxiliary variable that checks how many steps
    // left are there
    // int n_steps_left;
    // bool queuing_started   = false;

public:

    void operator+=(Agent<TSeq> * p);
    void operator-=(Agent<TSeq> * p);
    epiworld_fast_int & operator[](epiworld_fast_uint i);

    // void initialize(Model<TSeq> * m, Agent<TSeq> * p);
    void reset();

    bool operator==(const Queue<TSeq> & other) const;
    bool operator!=(const Queue<TSeq> & other) const {return !operator==(other);};

    static const int NoOne    = 0;
    static const int OnlySelf = 1;
    static const int Everyone = 2;

};

template<typename TSeq>
inline void Queue<TSeq>::operator+=(Agent<TSeq> * p)
{

    if (++active[p->id] == 1)
        n_in_queue++;

    for (auto n : (*p->neighbors))
    {

        if (++active[n] == 1)
            n_in_queue++;

    }

}

template<typename TSeq>
inline void Queue<TSeq>::operator-=(Agent<TSeq> * p)
{

    if (--active[p->id] == 0)
        n_in_queue--;

    for (auto n : (*p->neighbors))
    {
        if (--active[n] == 0)
            n_in_queue--;
    }

}

template<typename TSeq>
inline epiworld_fast_int & Queue<TSeq>::operator[](epiworld_fast_uint i)
{
    return active[i];
}

template<typename TSeq>
inline void Queue<TSeq>::reset()
{

    if (n_in_queue)
    {

        for (auto & q : this->active)
            q = 0;

        n_in_queue = 0;
        
    }

    active.resize(model->size(), 0);

}

template<typename TSeq>
inline bool Queue<TSeq>::operator==(const Queue<TSeq> & other) const 
{
    if (active.size() != other.active.size())
        return false;

    for (size_t i = 0u; i < active.size(); ++i)
    {
        if (active[i] != other.active[i])
            return false;
    }

    return true;
}

#endif