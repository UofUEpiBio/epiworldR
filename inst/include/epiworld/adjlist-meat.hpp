#ifndef EPIWORLD_ADJLIST_MEAT_HPP
#define EPIWORLD_ADJLIST_MEAT_HPP

inline AdjList::AdjList(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) : directed(directed) {


    dat.resize(size, std::map<int,int>({}));
    int max_id = size - 1;

    int i,j;
    for (int m = 0; m < static_cast<int>(source.size()); ++m)
    {

        i = source[m];
        j = target[m];

        if (i > max_id)
            throw std::range_error(
                "The source["+std::to_string(m)+"] = " + std::to_string(i) +
                " is above the max_id " + std::to_string(max_id)
                );

        if (j > max_id)
            throw std::range_error(
                "The target["+std::to_string(m)+"] = " + std::to_string(j) +
                " is above the max_id " + std::to_string(max_id)
                );

        // Adding nodes
        if (dat[i].find(j) == dat[i].end())
            dat[i].insert(std::pair<int, int>(j, 1u));
        else
            dat[i][j]++; 
        
        if (!directed)
        {

            if (dat[j].find(i) == dat[j].end())
                dat[j].insert(std::pair<int, int>(i, 1u));
            else
                dat[j][i]++;

        }

        E++;

    }

    N = size;

    return;

}


inline AdjList::AdjList(AdjList && a) :
    dat(std::move(a.dat)),
    directed(a.directed),
    N(a.N),
    E(a.E)
{

}

inline AdjList::AdjList(const AdjList & a) :
    dat(a.dat),
    directed(a.directed),
    N(a.N),
    E(a.E)
{

}

inline AdjList& AdjList::operator=(const AdjList& a)
{
    if (this == &a)
        return *this;

    this->dat = a.dat;
    this->directed = a.directed;
    this->N = a.N;
    this->E = a.E;

    return *this;
}

inline void AdjList::read_edgelist(
    std::string fn,
    int size,
    int skip,
    bool directed
) {

    int i,j;
    std::ifstream filei(fn);

    if (!filei)
        throw std::logic_error("The file " + fn + " was not found.");

    int linenum = 0;
    std::vector< int > source_;
    std::vector< int > target_;

    source_.reserve(1e5);
    target_.reserve(1e5);

    int max_id = size - 1;

    while (!filei.eof())
    {

        if (linenum++ < skip)
            continue;

        filei >> i >> j;

        // Looking for exceptions
        if (filei.bad())
            throw std::logic_error(
                "I/O error while reading the file " +
                fn
            );

        if (filei.fail())
            break;

        if (i > max_id)
            throw std::range_error(
                "The source["+std::to_string(linenum)+"] = " + std::to_string(i) +
                " is above the max_id " + std::to_string(max_id)
                );

        if (j > max_id)
            throw std::range_error(
                "The target["+std::to_string(linenum)+"] = " + std::to_string(j) +
                " is above the max_id " + std::to_string(max_id)
                );

        source_.push_back(i);
        target_.push_back(j);

    }

    // Now using the right constructor
    *this = AdjList(source_, target_, size, directed);

    return;

}

inline std::map<int,int> AdjList::operator()(
    epiworld_fast_uint i
    ) const {

    if (i >= N)
        throw std::range_error(
            "The vertex id " + std::to_string(i) + " is not in the network."
            );

    return dat[i];

}

inline void AdjList::print(epiworld_fast_uint limit) const {


    epiworld_fast_uint counter = 0;
    printf_epiworld("Nodeset:\n");
    int i = -1;
    for (auto & n : dat)
    {

        if (counter++ > limit)
            break;

        printf_epiworld("  % 3i: {", ++i);
        int niter = 0;
        for (auto n_n : n)
            if (++niter < static_cast<int>(n.size()))
            {    
                printf_epiworld("%i, ", static_cast<int>(n_n.first));
            }
            else {
                printf_epiworld("%i}\n", static_cast<int>(n_n.first));
            }
    }

    if (limit < dat.size())
    {
        printf_epiworld(
            "  (... skipping %i records ...)\n",
            static_cast<int>(dat.size() - limit)
            );
    }

}

inline size_t AdjList::vcount() const 
{
    return N;
}

inline size_t AdjList::ecount() const 
{
    return E;
}

inline bool AdjList::is_directed() const {

    if (dat.size() == 0u)
        throw std::logic_error("The edgelist is empty.");
    
    return directed;
    
}

#endif