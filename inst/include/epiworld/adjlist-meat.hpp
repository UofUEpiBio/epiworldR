#ifndef EPIWORLD_ADJLIST_MEAT_HPP
#define EPIWORLD_ADJLIST_MEAT_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <stdexcept>
#include <fstream>
#include "config.hpp"
#include "adjlist-bones.hpp"

inline AdjList::AdjList(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) : directed(directed) {

    if (size < 0)
        throw std::length_error(
            "The size of the network cannot be negative (" +
            std::to_string(size) + ")."
            );

    if (source.size() != target.size())
        throw std::length_error(
            "source and target must have the same length (" +
            std::to_string(source.size()) + " vs " +
            std::to_string(target.size()) + ")."
            );

    N = size;
    E = source.size();
    int max_id = size - 1;

    // Counting sort of the tie ends by vertex: first the size of each row (an
    // undirected tie counts in both), then the rows themselves
    start.assign(N + 1u, 0u);
    for (size_t m = 0u; m < source.size(); ++m)
    {

        int i = source[m];
        int j = target[m];

        if ((i < 0) || (i > max_id))
            throw std::range_error(
                "The source["+std::to_string(m)+"] = " + std::to_string(i) +
                " is above the max_id " + std::to_string(max_id)
                );

        if ((j < 0) || (j > max_id))
            throw std::range_error(
                "The target["+std::to_string(m)+"] = " + std::to_string(j) +
                " is above the max_id " + std::to_string(max_id)
                );

        start[i + 1]++;
        if (!directed)
            start[j + 1]++;

    }

    std::partial_sum(start.begin(), start.end(), start.begin());

    ids.resize(start[N]);
    std::vector< size_t > next(start.begin(), start.end() - 1);
    for (size_t m = 0u; m < source.size(); ++m)
    {
        ids[next[source[m]]++] = target[m];
        if (!directed)
            ids[next[target[m]]++] = source[m];
    }

    // Each row sorted and without repeats, moved down over the repeats
    // dropped before it
    size_t kept = 0u;
    for (size_t i = 0u; i < N; ++i)
    {

        auto first = ids.begin() + static_cast< std::ptrdiff_t >(start[i]);
        auto last  = ids.begin() + static_cast< std::ptrdiff_t >(start[i + 1u]);

        std::sort(first, last);
        last = std::unique(first, last);

        start[i] = kept;
        for (auto it = first; it != last; ++it)
            ids[kept++] = *it;

    }

    start[N] = kept;
    ids.resize(kept);

    return;

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

inline std::vector< int > AdjList::operator()(
    epiworld_fast_uint i
    ) const {

    if (i >= N)
        throw std::range_error(
            "The vertex id " + std::to_string(i) + " is not in the network."
            );

    return std::vector< int >(
        ids.begin() + static_cast< std::ptrdiff_t >(start[i]),
        ids.begin() + static_cast< std::ptrdiff_t >(start[i + 1u])
    );

}

inline void AdjList::print(epiworld_fast_uint limit) const {

    printf_epiworld("Nodeset:\n");
    for (epiworld_fast_uint i = 0u; (i < N) && (i <= limit); ++i)
    {

        printf_epiworld("  % 3i: {", static_cast<int>(i));
        for (size_t k = start[i]; k < start[i + 1u]; ++k)
        {
            printf_epiworld("%s%i", (k > start[i]) ? ", " : "", ids[k]);
        }
        printf_epiworld("}\n");

    }

    if (limit < N)
    {
        printf_epiworld(
            "  (... skipping %i records ...)\n",
            static_cast<int>(N - limit)
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

    if (N == 0u)
        throw std::logic_error("The edgelist is empty.");
    
    return directed;
    
}

#endif