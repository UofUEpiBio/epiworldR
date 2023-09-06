#ifndef EPIWORLD_NETWORK_BONES_HPP
#define EPIWORLD_NETWORK_BONES_HPP
template<typename Nettype, typename Nodetype, typename Edgetype>
class Network {

private:
    NType data;

public:

    NType();

    Edgetype operator()(int i, int j);

    bool is_directed() const;
    size_t vcount() const;
    size_t ecount() const;

    void add_edge(int i, int j);
    void rm_edge(int i, int j);
    

};

#endif