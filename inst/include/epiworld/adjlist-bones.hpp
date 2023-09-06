#ifndef EPIWORLD_ADJLIST_BONES_HPP
#define EPIWORLD_ADJLIST_BONES_HPP

class AdjList {
private:

    std::vector<std::map<int, int>> dat;
    bool directed;
    epiworld_fast_uint N = 0;
    epiworld_fast_uint E = 0;

public:

    AdjList() {};

    /**
     * @brief Construct a new Adj List object
     * 
     * @details 
     * Ids in the network are assume to range from `0` to `size - 1`.
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

    AdjList(AdjList && a); // Move constructor
    AdjList(const AdjList & a); // Copy constructor
    AdjList& operator=(const AdjList& a);


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

    std::map<int, int> operator()(
        epiworld_fast_uint i
        ) const;
        
    void print(epiworld_fast_uint limit = 20u) const;
    size_t vcount() const; ///< Number of vertices/nodes in the network.
    size_t ecount() const; ///< Number of edges/arcs/ties in the network.
    
    std::vector<std::map<int,int>> & get_dat() {
        return dat;
    };

    bool is_directed() const; ///< `true` if the network is directed.

};


#endif

