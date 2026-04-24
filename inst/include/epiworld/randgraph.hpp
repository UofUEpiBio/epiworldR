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
        _degree0[i] = model->get_agents()[i].get_neighbors(*model).size();
    #endif

    // Identifying individuals with degree > 0
    std::vector< epiworld_fast_uint > non_isolates;
    std::vector< epiworld_double > weights;
    epiworld_double nedges = 0.0;

    for (epiworld_fast_uint i = 0u; i < agents->size(); ++i)
    {
        if (agents->operator[](i).get_neighbors(*model).size() > 0u)
        {
            non_isolates.push_back(i);
            epiworld_double wtemp = static_cast<epiworld_double>(
                agents->operator[](i).get_neighbors(*model).size()
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
        int id01 = model->runif_index(p0.get_n_neighbors());
        int id11 = model->runif_index(p1.get_n_neighbors());

        // Get the actual neighbor IDs that will be swapped
        auto neighbors_p0 = p0.get_neighbors(*model);
        auto neighbors_p1 = p1.get_neighbors(*model);
        size_t neighbor_id_01 = neighbors_p0[id01]->get_id();
        size_t neighbor_id_11 = neighbors_p1[id11]->get_id();

        // Check if the swap would create self-loops or invalid configurations
        // After swap: p0 will be connected to neighbor_id_11, p1 to neighbor_id_01
        // Skip if:
        // 1. neighbor_id_01 == neighbor_id_11 (swapping the same neighbor)
        // 2. neighbor_id_01 == non_isolates[id1] (p1's new neighbor would be p1 itself)
        // 3. neighbor_id_11 == non_isolates[id0] (p0's new neighbor would be p0 itself)
        if (neighbor_id_01 == neighbor_id_11 ||
            neighbor_id_01 == non_isolates[id1] ||
            neighbor_id_11 == non_isolates[id0]) {
            continue;
        }

        // Check if the swap would create duplicate edges
        // After swap: p0 will be connected to neighbor_id_11, p1 to neighbor_id_01
        bool would_create_duplicate = false;
        for (auto* n : neighbors_p0) {
            if (n->get_id() == neighbor_id_11 && n != neighbors_p0[id01]) {
                would_create_duplicate = true;
                break;
            }
        }
        if (!would_create_duplicate) {
            for (auto* n : neighbors_p1) {
                if (n->get_id() == neighbor_id_01 && n != neighbors_p1[id11]) {
                    would_create_duplicate = true;
                    break;
                }
            }
        }

        if (would_create_duplicate) {
            continue; // Skip this rewire attempt
        }

        // When rewiring, we need to flip the individuals from the other
        // end as well, since we are dealing withi an undirected graph
        
        // Swap neighbors between the two agents
        // Note: id0 and id1 are indices in non_isolates, not agent IDs
        model->get_agents()[non_isolates[id0]].swap_neighbors(
            model->get_agents()[non_isolates[id1]],
            id01,
            id11,
            *model
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
    
    std::vector< int > non_isolates;
    non_isolates.reserve(nties.size());

    std::vector< epiworld_double > weights;
    weights.reserve(nties.size());

    epiworld_double nedges = 0.0;
    auto & dat = agents->get_dat();

    if (dat.size() > nties.size())
        throw std::logic_error("Inconsistent adjacency list data.");

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
        int id01 = model->runif_index(p0.size());
        int id11 = model->runif_index(p1.size());

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

        // When rewiring, we need to actually swap the edges, not just the weights
        // We'll swap edges: (id0, id01) <-> (id1, id11)
        // After swap: (id0, id11) and (id1, id01)
        // But first, check if the swap would create duplicate or self-loop edges
        
        // Check for self-loops (new edge would connect node to itself)
        if (id01 == non_isolates[id1] || id11 == non_isolates[id0]) {
            continue;
        }
        
        // Check for duplicate edges (new edge already exists)
        if (p0.find(id11) != p0.end() || p1.find(id01) != p1.end()) {
            continue; // Skip this rewire attempt to avoid duplicate edges
        }
        
        // Check if we're trying to swap the same neighbor
        if (id01 == id11) {
            continue;
        }
        
        // Save the weights before removing edges
        int weight_0_01 = p0[id01];
        int weight_1_11 = p1[id11];
        
        // Remove old edges from ego perspectives
        p0.erase(id01);
        p1.erase(id11);
        
        // Add new edges from ego perspectives with swapped alters
        p0[id11] = weight_0_01;
        p1[id01] = weight_1_11;
        
        // For undirected graphs, also update from alter perspectives
        if (!agents->is_directed())
        {
            std::map<int,int> & p01 = agents->get_dat()[id01];
            std::map<int,int> & p11 = agents->get_dat()[id11];
            
            // Save weights from alter perspectives
            int weight_01_0 = p01[non_isolates[id0]];
            int weight_11_1 = p11[non_isolates[id1]];
            
            // Remove old edges from alter perspectives
            p01.erase(non_isolates[id0]);
            p11.erase(non_isolates[id1]);
            
            // Add new edges from alter perspectives
            p01[non_isolates[id1]] = weight_11_1;
            p11[non_isolates[id0]] = weight_01_0;
        }

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

/**
 * Threshold for the hybrid row-advance in Batagelj-Brandes.
 *
 * When a geometric skip lands past the current row of the lower-triangle
 * enumeration, the flat index must be mapped back to (i, j). For small
 * excess (j - i <= threshold) the integer while-loop is fastest; for
 * large excess the O(1) inverse-triangular-number formula avoids
 * O(sqrt(excess)) iterations. 200 keeps the loop body to at most ~20
 * iterations (sqrt(2 * 200) ~ 20) before switching to FPU math.
 */
static constexpr long long BB_FPU_THRESHOLD = 200;

/**
 * @brief Generates an Erdős–Rényi random graph G(n, p) using the
 * Batagelj-Brandes algorithm.
 *
 * Uses geometric skipping to advance through the space of possible edges,
 * generating each edge independently with probability p. This produces
 * **unbiased** edge counts in O(n + m) expected time, where m is the
 * number of edges generated.
 *
 * Reference: V. Batagelj and U. Brandes, "Efficient generation of large
 * random networks," Physical Review E, vol. 71, no. 3, 2005.
 *
 * @tparam TSeq Type of the sequence (template parameter of the model).
 * @param n Number of nodes.
 * @param p Edge probability.
 * @param directed Whether the graph is directed.
 * @param model A reference to the Model, used for random number generation.
 * @return An AdjList representing the generated network.
 */
template<typename TSeq>
inline AdjList rgraph_bernoulli(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    std::vector< int > source;
    std::vector< int > target;

    // Pre-allocate based on expected edge count
    double n_possible = directed
        ? static_cast<double>(n) * (n - 1.0)
        : static_cast<double>(n) * (n - 1.0) / 2.0;
    source.reserve(static_cast<size_t>(n_possible * p * 1.1));
    target.reserve(static_cast<size_t>(n_possible * p * 1.1));

    // Batagelj-Brandes: use geometric skipping to enumerate edges.
    // log(1 - p) is precomputed; each skip is 1 + floor(log(U) / lp).
    double lp = std::log(1.0 - p);

    if (!directed)
    {
        // Undirected: iterate over (i, j) with i > j
        long long i = 1;
        long long j = -1;

        while (i < static_cast<long long>(n))
        {
            // Geometric skip
            double r = model.runif();

            // Avoid log(0)
            if (r == 0.0)
                r = std::numeric_limits<double>::min();

            j += 1 + static_cast<long long>(
                std::floor(std::log(r) / lp)
            );

            // Advance to the next row if j >= i.
            // Hybrid approach: for large skips, use O(1) FPU
            // inverse-triangular-number formula instead of the
            // integer loop, avoiding O(sqrt(skip)) iterations.
            if (j >= i && i < static_cast<long long>(n))
            {
                if (j - i > BB_FPU_THRESHOLD)
                {
                    // O(1) FPU mapping: flat index → (i, j) via
                    // inverse triangular number. `flat` fits in
                    // long long since n ≤ INT_MAX (AdjList limit),
                    // so i*(i-1)/2 + j < n*(n-1)/2 ≪ 2^62.
                    long long flat = i * (i - 1) / 2 + j;
                    double df = static_cast<double>(flat);
                    i = static_cast<long long>(std::floor(
                        (1.0 + std::sqrt(1.0 + 8.0 * df)) / 2.0
                    ));
                    j = flat - i * (i - 1) / 2;
                    // Floating-point correction: sqrt may round i
                    // up or down by 1 for large flat values, so j
                    // may land outside [0, i). Each correction loop
                    // runs at most 1–2 iterations; only one of the
                    // two can execute for a given (i, j).
                    while (j >= i) { j -= i; i++; }
                    while (j < 0)  { i--; j += i; }
                }
                else
                {
                    while (j >= i && i < static_cast<long long>(n))
                    {
                        j -= i;
                        i++;
                    }
                }
            }

            if (i < static_cast<long long>(n))
            {
                source.push_back(static_cast<int>(i));
                target.push_back(static_cast<int>(j));
            }
        }
    }
    else
    {
        // Directed: iterate over (i, j) with i != j
        // We linearize pairs: for row i, columns [0..n-1] excluding i.
        // Position k in [0, n*(n-1)) maps to row = k/(n-1),
        // col = k%(n-1), if col >= row then col++.
        long long total = static_cast<long long>(n) * (n - 1);
        long long k = -1;

        while (true)
        {
            double r = model.runif();
            if (r == 0.0)
                r = std::numeric_limits<double>::min();

            k += 1 + static_cast<long long>(
                std::floor(std::log(r) / lp)
            );

            if (k >= total)
                break;

            long long row = k / static_cast<long long>(n - 1);
            long long col = k % static_cast<long long>(n - 1);
            if (col >= row)
                col++;

            source.push_back(static_cast<int>(row));
            target.push_back(static_cast<int>(col));
        }
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

/**
 * @brief Generates a network using a Stochastic Block Model (SBM).
 *
 * This function creates a random undirected network where connections between
 * agents are determined by a mixing matrix and block (group) membership.
 * Each potential edge between agents in blocks \f$g\f$ and \f$h\f$ is included
 * independently with probability \f$M(g, h) / n_h\f$, where \f$M(g, h)\f$ is
 * the entry of the mixing matrix and \f$n_h\f$ is the number of agents in
 * block \f$h\f$.
 *
 * The mixing matrix is **not** row-stochastic. For undirected networks, its
 * row sums should be interpreted as the average expected degree for agents in
 * that group only when the balance condition
 * \f$M(g, h) \times n_g = M(h, g) \times n_h\f$ holds for all pairs
 * \f$(g, h)\f$. In that balanced case,
 * \f[
 *   \mathbb{E}[\text{degree of agent in group } g] = \sum_{h} M(g, h)
 * \f]
 *
 * For undirected networks, the implementation samples between-block edges only
 * once for each unordered pair of blocks (with \f$g \le h\f$) using
 * \f$M(g, h) / n_h\f$. Therefore, if the balance condition is violated, the
 * row sums generally do not match the expected degree from every group's
 * perspective.
 *
 * **Unbiased generation.** Edge generation uses the Batagelj-Brandes
 * algorithm (geometric skipping) independently for each block pair.
 * Each possible edge is included with exactly the correct probability,
 * with no duplicate edge collisions. The expected number of edges for a
 * block pair is exactly:
 * - Within-block: \f$\binom{n_g}{2} \times p_{gg}\f$
 * - Between-block: \f$n_g \times n_h \times p_{gh}\f$
 *
 * **Balance condition warning.** If the balance condition is violated for
 * any pair of blocks and the model's verbose mode is on, a warning is
 * printed indicating which entries are being silently ignored.
 *
 * Reference: V. Batagelj and U. Brandes, "Efficient generation of large
 * random networks," Physical Review E, vol. 71, no. 3, 2005.
 *
 * @tparam TSeq Type of the sequence (template parameter of the model).
 * @param block_sizes A vector of size \f$K\f$ indicating the number of agents
 *   per block. The total number of agents is \f$\sum_k \text{block\_sizes}[k]\f$.
 * @param mixing_matrix A vector of size \f$K \times K\f$ representing the
 *   mixing matrix. The entry \f$M(g, h)\f$ controls the expected number of
 *   connections agents in group \f$g\f$ make with agents in group \f$h\f$.
 * @param row_major If `true`, the mixing matrix is stored in row-major order
 *   (C-style), i.e., \f$M(g, h) = \text{mixing\_matrix}[g \times K + h]\f$.
 *   If `false`, column-major order (Fortran-style) is assumed, i.e.,
 *   \f$M(g, h) = \text{mixing\_matrix}[h \times K + g]\f$.
 * @param model A reference to the Model, used for random number generation.
 * @return An AdjList representing the generated undirected network.
 *
 * @throws std::length_error If `block_sizes` is empty, `mixing_matrix` size
 *   does not equal \f$K^2\f$, or any block size is zero.
 * @throws std::range_error If any mixing matrix entry is negative or if
 *   \f$M(g, h) / n_h > 1\f$ for any pair of blocks.
 */
template<typename TSeq>
inline AdjList rgraph_sbm(
    const std::vector< size_t > & block_sizes,
    const std::vector< double > & mixing_matrix,
    bool row_major,
    Model<TSeq> & model
) {

    size_t n_blocks = block_sizes.size();

    if (n_blocks == 0u)
        throw std::length_error(
            "block_sizes must have at least one element."
        );

    if (mixing_matrix.size() != n_blocks * n_blocks)
        throw std::length_error(
            "mixing_matrix must have size K * K, where K = block_sizes.size(). "
            "Expected " + std::to_string(n_blocks * n_blocks) +
            " but got " + std::to_string(mixing_matrix.size()) + "."
        );

    // Total number of agents
    size_t n = 0u;
    const size_t max_adjlist_n = static_cast<size_t>(
        std::numeric_limits<int>::max()
    );
    for (auto bs : block_sizes)
    {
        if (bs == 0u)
            throw std::length_error(
                "All block sizes must be positive."
            );

        if (bs > max_adjlist_n - n)
            throw std::length_error(
                "The total number of agents implied by block_sizes exceeds "
                "the maximum supported graph size (" +
                std::to_string(std::numeric_limits<int>::max()) + ")."
            );

        n += bs;
    }

    // Compute the starting index for each block
    std::vector< size_t > block_start(n_blocks, 0u);
    for (size_t g = 1u; g < n_blocks; ++g)
        block_start[g] = block_start[g - 1u] + block_sizes[g - 1u];

    // Validate mixing matrix entries and compute probabilities.
    // p(g, h) = M(g, h) / n_h for the undirected case.
    // We validate that all entries are non-negative and p(g,h) <= 1.
    // Also estimate the total expected edge count for vector reservation.
    double estimated_edges = 0.0;
    for (size_t g = 0u; g < n_blocks; ++g)
    {
        for (size_t h = 0u; h < n_blocks; ++h)
        {
            double m_gh = row_major
                ? mixing_matrix[g * n_blocks + h]
                : mixing_matrix[h * n_blocks + g];

            if (m_gh < 0.0)
                throw std::range_error(
                    "Mixing matrix entries must be non-negative. "
                    "Entry (" + std::to_string(g) + ", " +
                    std::to_string(h) + ") = " + std::to_string(m_gh) + "."
                );

            double p_gh = m_gh / static_cast<double>(block_sizes[h]);
            if (p_gh > 1.0)
                throw std::range_error(
                    "Mixing matrix entry (" + std::to_string(g) + ", " +
                    std::to_string(h) + ") = " + std::to_string(m_gh) +
                    " implies probability " + std::to_string(p_gh) +
                    " > 1 (block size = " +
                    std::to_string(block_sizes[h]) + ")."
                );
        }

        // Accumulate expected edges for upper triangle (g <= h)
        for (size_t h = g; h < n_blocks; ++h)
        {
            double m_gh = row_major
                ? mixing_matrix[g * n_blocks + h]
                : mixing_matrix[h * n_blocks + g];

            double p_gh = m_gh / static_cast<double>(block_sizes[h]);
            if (p_gh <= 0.0)
                continue;

            double n_possible;
            if (g == h)
                n_possible = static_cast<double>(block_sizes[g]) *
                    static_cast<double>(block_sizes[g] - 1u) / 2.0;
            else
                n_possible = static_cast<double>(block_sizes[g]) *
                    static_cast<double>(block_sizes[h]);

            estimated_edges += n_possible * std::min(p_gh, 1.0);
        }
    }

    // Warn when the balance condition M(g,h)*n_g = M(h,g)*n_h is
    // violated for any pair. The upper triangle is used for sampling,
    // so the lower-triangle entries are silently ignored.
    if (model.get_verbose())
    {
        for (size_t g = 0u; g < n_blocks; ++g)
        {
            for (size_t h = g + 1u; h < n_blocks; ++h)
            {
                double m_gh = row_major
                    ? mixing_matrix[g * n_blocks + h]
                    : mixing_matrix[h * n_blocks + g];
                double m_hg = row_major
                    ? mixing_matrix[h * n_blocks + g]
                    : mixing_matrix[g * n_blocks + h];

                double lhs = m_gh * static_cast<double>(block_sizes[g]);
                double rhs = m_hg * static_cast<double>(block_sizes[h]);

                if (std::abs(lhs - rhs) > 1e-10 *
                    std::max(std::abs(lhs), std::abs(rhs)))
                {
                    printf_epiworld(
                        "Warning [rgraph_sbm]: balance condition violated "
                        "for blocks (%zu, %zu): M(%zu,%zu)*n_%zu = %.4f "
                        "!= M(%zu,%zu)*n_%zu = %.4f. Only M(%zu,%zu)/n_%zu "
                        "is used for edge sampling.\n",
                        g, h,
                        g, h, g, lhs,
                        h, g, h, rhs,
                        g, h, h
                    );
                }
            }
        }
    }

    std::vector< int > source;
    std::vector< int > target;

    // Pre-allocate based on estimated edge count to avoid
    // repeated reallocations during push_back.
    source.reserve(static_cast<size_t>(estimated_edges * 1.1));
    target.reserve(static_cast<size_t>(estimated_edges * 1.1));

    // Batagelj-Brandes geometric skipping for each block pair.
    // For each pair (g, h) with g <= h, we enumerate possible edges
    // and skip ahead using geometric random variables, producing each
    // edge independently with exactly probability p_gh. No duplicates,
    // no bias.
    for (size_t g = 0u; g < n_blocks; ++g)
    {
        for (size_t h = g; h < n_blocks; ++h)
        {

            double m_gh = row_major
                ? mixing_matrix[g * n_blocks + h]
                : mixing_matrix[h * n_blocks + g];

            // Probability of an edge between agents in blocks g and h
            double p_gh = m_gh / static_cast<double>(block_sizes[h]);

            if (p_gh <= 0.0)
                continue;

            // Fast path: p_gh >= 1.0 means all edges are present
            if (p_gh >= 1.0)
            {
                if (g == h)
                {
                    // All unique pairs within the block
                    for (size_t i = block_start[g];
                         i < block_start[g] + block_sizes[g]; ++i)
                    {
                        for (size_t j = i + 1u;
                             j < block_start[g] + block_sizes[g]; ++j)
                        {
                            source.push_back(static_cast<int>(i));
                            target.push_back(static_cast<int>(j));
                        }
                    }
                }
                else
                {
                    // All pairs between blocks g and h
                    for (size_t i = block_start[g];
                         i < block_start[g] + block_sizes[g]; ++i)
                    {
                        for (size_t j = block_start[h];
                             j < block_start[h] + block_sizes[h]; ++j)
                        {
                            source.push_back(static_cast<int>(i));
                            target.push_back(static_cast<int>(j));
                        }
                    }
                }
                continue;
            }

            double lp = std::log(1.0 - p_gh);

            if (g == h)
            {
                // Within-block: iterate (i, j) with i > j, both in
                // [0, n_g). Use BB geometric skipping.
                long long n_g = static_cast<long long>(block_sizes[g]);
                long long i = 1;
                long long j = -1;

                while (i < n_g)
                {
                    double r = model.runif();
                    if (r == 0.0)
                        r = std::numeric_limits<double>::min();

                    j += 1 + static_cast<long long>(
                        std::floor(std::log(r) / lp)
                    );

                    // Hybrid row-advance (see BB_FPU_THRESHOLD
                    // comment in rgraph_bernoulli above).
                    if (j >= i && i < n_g)
                    {
                        if (j - i > BB_FPU_THRESHOLD)
                        {
                            long long flat = i * (i - 1) / 2 + j;
                            double df = static_cast<double>(flat);
                            i = static_cast<long long>(std::floor(
                                (1.0 + std::sqrt(1.0 + 8.0 * df)) / 2.0
                            ));
                            j = flat - i * (i - 1) / 2;
                            while (j >= i) { j -= i; i++; }
                            while (j < 0)  { i--; j += i; }
                        }
                        else
                        {
                            while (j >= i && i < n_g)
                            {
                                j -= i;
                                i++;
                            }
                        }
                    }

                    if (i < n_g)
                    {
                        source.push_back(
                            static_cast<int>(block_start[g] +
                                static_cast<size_t>(i))
                        );
                        target.push_back(
                            static_cast<int>(block_start[g] +
                                static_cast<size_t>(j))
                        );
                    }
                }
            }
            else
            {
                // Between-block: iterate (i, j) with i in [0, n_g),
                // j in [0, n_h). Linearize as k = i * n_h + j.
                long long n_g = static_cast<long long>(block_sizes[g]);
                long long n_h = static_cast<long long>(block_sizes[h]);
                long long total = n_g * n_h;
                long long k = -1;

                while (true)
                {
                    double r = model.runif();
                    if (r == 0.0)
                        r = std::numeric_limits<double>::min();

                    k += 1 + static_cast<long long>(
                        std::floor(std::log(r) / lp)
                    );

                    if (k >= total)
                        break;

                    long long row = k / n_h;
                    long long col = k % n_h;

                    source.push_back(
                        static_cast<int>(block_start[g] +
                            static_cast<size_t>(row))
                    );
                    target.push_back(
                        static_cast<int>(block_start[h] +
                            static_cast<size_t>(col))
                    );
                }
            }

        }
    }

    return AdjList(source, target, static_cast<int>(n), false);

}

#endif
