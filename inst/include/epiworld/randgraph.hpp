#ifndef EPIWORLD_RANDGRA
#define EPIWORLD_RANDGRA

template<typename TSeq>
class Model;

template<typename TSeq>
class Agent;

class AdjList;


template<typename TSeq, typename TDat>
inline void rewire_degseq(
    TDat * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    );

template<typename TSeq = EPI_DEFAULT_TSEQ>
inline void rewire_degseq(
    std::vector< Agent<TSeq> > * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    )
{

    #ifdef EPI_DEBUG
    std::vector< int > _degree0(agents->size(), 0);
    for (size_t i = 0u; i < _degree0.size(); ++i)
        _degree0[i] = model->get_agents()[i].get_neighbors().size();
    #endif

    // Identifying individuals with degree > 0
    std::vector< epiworld_fast_uint > non_isolates;
    std::vector< epiworld_double > weights;
    epiworld_double nedges = 0.0;
    
    for (epiworld_fast_uint i = 0u; i < agents->size(); ++i)
    {
        if (agents->operator[](i).get_neighbors().size() > 0u)
        {
            non_isolates.push_back(i);
            epiworld_double wtemp = static_cast<epiworld_double>(
                agents->operator[](i).get_neighbors().size()
                );
            weights.push_back(wtemp);
            nedges += wtemp;
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (epiworld_fast_uint i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    // Only swap if needed
    epiworld_fast_uint N = non_isolates.size();
    epiworld_double prob;
    int nrewires = floor(proportion * nedges);
    while (nrewires-- > 0)
    {

        // Picking egos
        prob = model->runif();
        int id0 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id0 = i;
                break;
            }

        prob = model->runif();
        int id1 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id1 = i;
                break;
            }

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= static_cast<int>(N))
            id1 = 0;

        Agent<TSeq> & p0 = agents->operator[](non_isolates[id0]);
        Agent<TSeq> & p1 = agents->operator[](non_isolates[id1]);

        // Picking alters (relative location in their lists)
        // In this case, these are uniformly distributed within the list
        int id01 = std::floor(p0.get_n_neighbors() * model->runif());
        int id11 = std::floor(p1.get_n_neighbors() * model->runif());

        // When rewiring, we need to flip the individuals from the other
        // end as well, since we are dealing withi an undirected graph
        
        // Finding what neighbour is id0
        model->get_agents()[id0].swap_neighbors(
            model->get_agents()[id1],
            id01,
            id11
            );
        

    }

    #ifdef EPI_DEBUG
    for (size_t _i = 0u; _i < _degree0.size(); ++_i)
    {
        if (_degree0[_i] != static_cast<int>(model->get_agents()[_i].get_n_neighbors()))
            throw std::logic_error("[epi-debug] Degree does not match afted rewire_degseq.");
    }
    #endif

    return;

}

template<typename TSeq>
inline void rewire_degseq(
    AdjList * agents,
    Model<TSeq> * model,
    epiworld_double proportion
    )
{

    // Identifying individuals with degree > 0
    std::vector< epiworld_fast_int > nties(agents->vcount(), 0); 

    #ifdef EPI_DEBUG
    std::vector< int > _degree0(agents->vcount(), 0);
    for (size_t i = 0u; i < _degree0.size(); ++i)
        _degree0[i] = agents->get_dat()[i].size();
    #endif
    
    std::vector< epiworld_fast_uint > non_isolates;
    non_isolates.reserve(nties.size());

    std::vector< epiworld_double > weights;
    weights.reserve(nties.size());

    epiworld_double nedges = 0.0;
    auto & dat = agents->get_dat();

    for (size_t i = 0u; i < dat.size(); ++i)
        nties[i] += dat[i].size();
    
    bool directed = agents->is_directed();
    for (size_t i = 0u; i < dat.size(); ++i)
    {
        if (nties[i] > 0)
        {
            non_isolates.push_back(i);
            if (directed)
            {
                weights.push_back( 
                    static_cast<epiworld_double>(nties[i])
                );
                nedges += static_cast<epiworld_double>(nties[i]);
            }
            else {
                weights.push_back( 
                    static_cast<epiworld_double>(nties[i])/2.0
                );
                nedges += static_cast<epiworld_double>(nties[i]) / 2.0;
            }
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (epiworld_fast_uint i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    // Only swap if needed
    epiworld_fast_uint N = non_isolates.size();
    epiworld_double prob;
    int nrewires = floor(proportion * nedges / (
        agents->is_directed() ? 1.0 : 2.0
    ));

    while (nrewires-- > 0)
    {

        // Picking egos
        prob = model->runif();
        int id0 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id0 = i;
                break;
            }

        prob = model->runif();
        int id1 = N - 1;
        for (epiworld_fast_uint i = 0u; i < N; ++i)
            if (prob <= weights[i])
            {
                id1 = i;
                break;
            }

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= static_cast<int>(N))
            id1 = 0;

        std::map<int,int> & p0 = agents->get_dat()[non_isolates[id0]];
        std::map<int,int> & p1 = agents->get_dat()[non_isolates[id1]];

        // Picking alters (relative location in their lists)
        // In this case, these are uniformly distributed within the list
        int id01 = std::floor(p0.size() * model->runif());
        int id11 = std::floor(p1.size() * model->runif());

        // Since it is a map, we need to find the actual ids (positions)
        // are not good enough.
        int count = 0;
        for (auto & n : p0)
            if (count++ == id01)
                id01 = n.first;

        count = 0;
        for (auto & n : p1)
            if (count++ == id11)
                id11 = n.first;

        // When rewiring, we need to flip the individuals from the other
        // end as well, since we are dealing withi an undirected graph
        
        // Finding what neighbour is id0
        if (!agents->is_directed())
        {

            std::map<int,int> & p01 = agents->get_dat()[id01];
            std::map<int,int> & p11 = agents->get_dat()[id11];

            std::swap(p01[id0], p11[id1]);
            
        }

        // Moving alter first
        std::swap(p0[id01], p1[id11]);

    }

    #ifdef EPI_DEBUG
    for (size_t _i = 0u; _i < _degree0.size(); ++_i)
    {
        if (_degree0[_i] != static_cast<int>(agents->get_dat()[_i].size()))
            throw std::logic_error(
                "[epi-debug] Degree does not match afted rewire_degseq. " +
                std::string("Expected: ") + 
                std::to_string(_degree0[_i]) + 
                std::string(", observed: ") +
                std::to_string(agents->get_dat()[_i].size())
                );
    }
    #endif


    return;

}

template<typename TSeq>
inline AdjList rgraph_bernoulli(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    std::vector< int > source;
    std::vector< int > target;

    // Checking the density (how many)
    std::binomial_distribution<> d(
        n * (n - 1.0) / (directed ? 1.0 : 2.0),
        p
    );

    epiworld_fast_uint m = d(model.get_rand_endgine());

    source.resize(m);
    target.resize(m);

    epiworld_fast_uint a,b;
    for (epiworld_fast_uint i = 0u; i < m; ++i)
    {
        a = floor(model.runif() * n);

        if (!directed)
            b = floor(model.runif() * a);
        else
        {
            b = floor(model.runif() * n);
            if (b == a)
                b++;
            
            if (b >= n)
                b = 0u;
        }

        source[i] = static_cast<int>(a);
        target[i] = static_cast<int>(b);

    }

    AdjList al(source, target, static_cast<int>(n), directed);

    return al;
    
}

template<typename TSeq>
inline AdjList rgraph_bernoulli2(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    std::vector< int > source;
    std::vector< int > target;

    // Checking the density (how many)
    std::binomial_distribution<> d(
        n * (n - 1.0) / (directed ? 1.0 : 2.0),
        p
    );

    // Need to compensate for the possible number of diagonal
    // elements sampled. If n * n, then each diag element has
    // 1/(n^2) chance of sampling

    epiworld_fast_uint m = d(model.get_rand_endgine());

    source.resize(m);
    target.resize(m);

    double n2 = static_cast<double>(n * n);

    int loc,row,col;
    for (epiworld_fast_uint i = 0u; i < m; ++i)
    {
        loc = floor(model.runif() * n2);
        col = floor(static_cast<double>(loc)/static_cast<double>(n));
        row = loc - row * n;

        // Undirected needs to swap
        if (!directed && (col > row))
            std::swap(col, row);

        source[i] = row;
        target[i] = col;

    }

    AdjList al(source, target, static_cast<int>(n), directed);

    return al;
    
}

inline AdjList rgraph_ring_lattice(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    bool directed = false
) {

    if ((n - 1u) < k)
        throw std::logic_error("k can be at most n - 1.");

    std::vector< int > source;
    std::vector< int > target;

    if (!directed)
        if (k > 1u) k = static_cast< size_t >(floor(k / 2.0));

    for (size_t i = 0; i < n; ++i)
    {

        for (size_t j = 1u; j <= k; ++j)
        {

            // Next neighbor
            size_t l = i + j;
            if (l >= n) l = l - n;

            source.push_back(i);
            target.push_back(l);

        }

    }

    return AdjList(source, target, n, directed);

}

/**
 * @brief Smallworld network (Watts-Strogatz)
 * 
 * @tparam TSeq 
 * @param n 
 * @param k 
 * @param p 
 * @param directed 
 * @param model 
 * @return AdjList 
 */
template<typename TSeq>
inline AdjList rgraph_smallworld(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    // Creating the ring lattice
    AdjList ring = rgraph_ring_lattice(n,k,directed);
    
    // Rewiring and returning
    if (k > 0u)
        rewire_degseq(&ring, &model, p);
        
    return ring;

}

/**
 * @brief Generates a blocked network
 * 
 * Since block sizes and number of connections between blocks are fixed,
 * this routine is fully deterministic.
 * 
 * @tparam TSeq 
 * @param n Size of the network
 * @param blocksize Size of the block.
 * @param ncons Number of connections between blocks
 * @param model A model
 * @return AdjList 
 */
template<typename TSeq>
inline AdjList rgraph_blocked(
    epiworld_fast_uint n,
    epiworld_fast_uint blocksize,
    epiworld_fast_uint ncons,
    Model<TSeq>&
) {

    std::vector< int > source_;
    std::vector< int > target_;

    size_t i = 0u;
    size_t cum_node_count = 0u;
    while (i < n)
    {

        for (size_t j = 0; j < blocksize; ++j)
        {

            for (size_t k = 0; k < j; ++k)
            {
                // No loops
                if (k == j)
                    continue;

                // Exists the loop in case there are no more 
                // nodes available
                if ((i + k) >= n)
                    break;

                source_.push_back(static_cast<int>(j + i));
                target_.push_back(static_cast<int>(k + i));
            }

            // No more nodes left to build connections
            if (++cum_node_count >= n)
                break;
            
        }

        // Connections between this and the previou sone
        if (i != 0)
        {

            size_t max_cons = std::min(ncons, n - cum_node_count);

            // Generating the connections
            for (size_t j = 0u; j < max_cons; ++j)
            {

                source_.push_back(static_cast<int>(i + j - blocksize));
                target_.push_back(static_cast<int>(i + j));

            }
        }

        i += blocksize;
        
    }
        
    return AdjList(source_, target_, n, false);

}

#endif
