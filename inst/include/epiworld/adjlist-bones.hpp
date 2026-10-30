#ifndef EPIWORLD_ADJLIST_BONES_HPP
#define EPIWORLD_ADJLIST_BONES_HPP


/**
 * @brief A network, as the ties of each vertex
 *
 * @details Row `i` holds the neighbors of vertex `i` (in a directed network,
 * the targets of its ties) in ascending order, without repeats; an undirected
 * tie is in the rows of both ends. The rows are stored back to back
 * (compressed rows): row `i` is `get_ids()` from `get_start()[i]` to
 * `get_start()[i + 1]`. A graph cannot be edited once built; the network layer
 * that loads it can (see `layer::Network`).
 */
class AdjList {
private:

    std::vector< size_t > start; ///< Where each row starts in `ids` (`N + 1`)
    std::vector< int > ids;      ///< The rows, back to back
    bool directed = false;
    epiworld_fast_uint N = 0;
    epiworld_fast_uint E = 0;

public:

    AdjList() {};

    /**
     * @brief Construct a new Adj List object
     *
     * @details
     * Ids in the network are assume to range from `0` to `size - 1`. A tie
     * given more than once is kept once.
     *
     * @param source Unsigned int vector with the source
     * @param target Unsigned int vector with the target
     * @param size Number of vertices in the network.
     * @param directed Bool true if the network is directed
     */
    AdjList(
        const std::vector< int > & source,
        const std::vector< int > & target,
        int size,
        bool directed
        );

    /**
     * @brief Read an edgelist
     *
     * Ids in the network are assume to range from `0` to `size - 1`.
     *
     * @param fn Path to the file
     * @param skip Number of lines to skip (e.g., 1 if there's a header)
     * @param directed `true` if the network is directed
     * @param size Number of vertices in the network.
     */
    void read_edgelist(
        std::string fn,
        int size,
        int skip = 0,
        bool directed = true
        );

    /// The neighbors of vertex `i` (a copy of its row).
    std::vector< int > operator()(
        epiworld_fast_uint i
        ) const;

    void print(epiworld_fast_uint limit = 20u) const;
    size_t vcount() const; ///< Number of vertices/nodes in the network.
    size_t ecount() const; ///< Number of edges given, repeats included.

    /// Where each row starts in `get_ids()`, with the end as the last entry.
    const std::vector< size_t > & get_start() const { return start; }

    /// The rows, back to back (see the class description).
    const std::vector< int > & get_ids() const { return ids; }

    bool is_directed() const; ///< `true` if the network is directed.

};


#endif
