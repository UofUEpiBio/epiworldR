#ifndef EPIWORLD_RANDOM_GRAPH_HPP
#define EPIWORLD_RANDOM_GRAPH_HPP

class RandGraph {

private:
    std::shared_ptr< std::mt19937 > engine;
    std::shared_ptr< std::uniform_real_distribution<> > unifd;
    int N = 0;
    bool initialized = false;

public:

    RandGraph(int N_) : N(N_) {};

    void init(int s);
    void set_rand_engine(std::mt19937 & e);
    epiworld_double runif();

};


inline void RandGraph::init(int s) {

    if (!engine)
        engine = std::make_shared< std::mt19937 >();

    engine->seed(s);

    if (!unifd)
        unifd = std::make_shared< std::uniform_real_distribution<> >(0, 1);

    initialized = true;


}

inline void RandGraph::set_rand_engine(std::mt19937 & e)
{

    engine = std::make_shared< std::mt19937 >(e);

}

inline epiworld_double RandGraph::runif() {

    if (!initialized)
        throw std::logic_error("The object has not been initialized");

    return (*unifd)(engine);

}

#endif