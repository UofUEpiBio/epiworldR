#ifndef EPIWORLD_RANDGRA
#define EPIWORLD_RANDGRA

template<typename TSeq>
class Model;

template<typename TSeq>
class Agent;

class AdjList;


/**
 * @brief The first `i` with `prob <= cumulative[i]`, or the last index.
 *
 * @details The draw of an agent with probability proportional to its weight,
 * given the cumulative probabilities (non-decreasing); the same index as a
 * scan from the start, found by bisection.
 */
inline size_t pick_cumulative(
    const std::vector< epiworld_double > & cumulative,
    epiworld_double prob
)
{
    auto it = std::lower_bound(cumulative.begin(), cumulative.end(), prob);
    return (it == cumulative.end()) ?
        cumulative.size() - 1u :
        static_cast< size_t >(it - cumulative.begin());
}

/**
 * @brief Degree-preserving rewiring of a graph under construction (used by
 * `rgraph_smallworld()`).
 *
 * @details `rows[i]` holds the neighbors of `i` in ascending order, without
 * repeats (in an undirected graph, each tie in the rows of both ends), and is
 * kept that way. Picks two agents with probability proportional to their
 * degree and one tie of each, uniformly, and swaps them, skipping swaps that
 * would create self-ties or repeated ties. A network layer is rewired with
 * `graph::rewire_degseq()` instead.
 */
template<typename TSeq>
inline void rewire_degseq_rows(
    std::vector< std::vector< int > > & rows,
    bool directed,
    Model<TSeq> & model,
    epiworld_double proportion
    )
{

    #ifdef EPI_DEBUG
    std::vector< size_t > _degree0(rows.size(), 0u);
    for (size_t i = 0u; i < rows.size(); ++i)
        _degree0[i] = rows[i].size();
    #endif

    std::vector< int > non_isolates;
    non_isolates.reserve(rows.size());

    std::vector< epiworld_double > weights;
    weights.reserve(rows.size());

    epiworld_double nedges = 0.0;
    for (size_t i = 0u; i < rows.size(); ++i)
    {
        if (rows[i].empty())
            continue;

        non_isolates.push_back(static_cast< int >(i));
        const epiworld_double nties = static_cast< epiworld_double >(rows[i].size());
        if (directed)
        {
            weights.push_back(nties);
            nedges += nties;
        }
        else
        {
            weights.push_back(nties / 2.0);
            nedges += nties / 2.0;
        }
    }

    if (non_isolates.size() == 0u)
        throw std::logic_error("The graph is completely disconnected.");

    // Cumulative probs
    weights[0u] /= nedges;
    for (size_t i = 1u; i < non_isolates.size(); ++i)
    {
         weights[i] /= nedges;
         weights[i] += weights[i - 1u];
    }

    const int N = static_cast< int >(non_isolates.size());
    int nrewires = floor(proportion * nedges / (directed ? 1.0 : 2.0));

    auto has = [](const std::vector< int > & row, int j) {
        return std::binary_search(row.begin(), row.end(), j);
    };

    auto insert = [](std::vector< int > & row, int j) {
        row.insert(std::lower_bound(row.begin(), row.end(), j), j);
    };

    auto erase = [](std::vector< int > & row, int j) {
        row.erase(std::lower_bound(row.begin(), row.end(), j));
    };

    while (nrewires-- > 0)
    {

        // Picking egos
        int id0 = static_cast< int >(pick_cumulative(weights, model.runif()));
        int id1 = static_cast< int >(pick_cumulative(weights, model.runif()));

        // Correcting for under or overflow.
        if (id1 == id0)
            id1++;

        if (id1 >= N)
            id1 = 0;

        const int ego0 = non_isolates[id0];
        const int ego1 = non_isolates[id1];
        std::vector< int > & p0 = rows[ego0];
        std::vector< int > & p1 = rows[ego1];

        // Picking alters (relative location in their lists), uniformly
        int id01 = model.runif_index(p0.size());
        int id11 = model.runif_index(p1.size());

        // From positions to ids. The scan keeps comparing the count with the
        // id it found (so a later position can replace it), as the version on
        // std::map did; kept so that the same seed gives the same graph.
        int count = 0;
        for (int j : p0)
            if (count++ == id01)
                id01 = j;

        count = 0;
        for (int j : p1)
            if (count++ == id11)
                id11 = j;

        // After the swap: (ego0, id11) and (ego1, id01). Skip self-ties,
        // repeated ties, and swaps of the same neighbor.
        if (id01 == ego1 || id11 == ego0)
            continue;

        if (has(p0, id11) || has(p1, id01))
            continue;

        if (id01 == id11)
            continue;

        erase(p0, id01);
        erase(p1, id11);
        insert(p0, id11);
        insert(p1, id01);

        // For undirected graphs, also update from alter perspectives
        if (!directed)
        {
            erase(rows[id01], ego0);
            erase(rows[id11], ego1);
            insert(rows[id01], ego1);
            insert(rows[id11], ego0);
        }

    }

    #ifdef EPI_DEBUG
    for (size_t _i = 0u; _i < _degree0.size(); ++_i)
    {
        if (_degree0[_i] != rows[_i].size())
            throw std::logic_error(
                "[epi-debug] Degree does not match afted rewire_degseq. " +
                std::string("Expected: ") +
                std::to_string(_degree0[_i]) +
                std::string(", observed: ") +
                std::to_string(rows[_i].size())
                );
    }
    #endif

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
 * @param p Edge probability, in [0, 1].
 * @param directed Whether the graph is directed.
 * @param model A reference to the Model, used for random number generation.
 * @return An AdjList representing the generated network.
 * @throws std::range_error if `p` is not in [0, 1].
 */
template<typename TSeq>
inline AdjList rgraph_bernoulli(
    epiworld_fast_uint n,
    epiworld_double p,
    bool directed,
    Model<TSeq> & model
) {

    // Written so that NaN fails too
    if (!((p >= 0.0) && (p <= 1.0)))
        throw std::range_error(
            "The edge probability must be in [0, 1]; got " +
            std::to_string(p) + "."
        );

    std::vector< int > source;
    std::vector< int > target;

    // No ties, and nothing to draw: with p = 0, log(1 - p) below is 0 and
    // the first skip would be infinite
    if (p == 0.0)
        return AdjList(source, target, static_cast<int>(n), directed);

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

    if (k == 0u)
        return ring;

    // Rewired as rows, which stay sorted, then loaded back
    const auto & start = ring.get_start();
    const auto & ids = ring.get_ids();
    std::vector< std::vector< int > > rows(n);
    for (size_t i = 0u; i < n; ++i)
        rows[i].assign(
            ids.begin() + static_cast< std::ptrdiff_t >(start[i]),
            ids.begin() + static_cast< std::ptrdiff_t >(start[i + 1u])
        );

    rewire_degseq_rows(rows, directed, model, p);

    // Each undirected tie once
    std::vector< int > source;
    std::vector< int > target;
    source.reserve(ids.size());
    target.reserve(ids.size());
    for (size_t i = 0u; i < n; ++i)
        for (int j : rows[i])
            if (directed || (static_cast< int >(i) <= j))
            {
                source.push_back(static_cast< int >(i));
                target.push_back(j);
            }

    return AdjList(source, target, static_cast< int >(n), directed);

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
