#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP agents_smallworld_cpp(
    SEXP m,
    unsigned int n = 1000,
    unsigned int k = 5,
    bool d = false,
    double p = .01

) {

  external_pointer<Model<>> ptr(m);
  ptr->agents_smallworld(n, k, d, p);

  return m;

}

[[cpp11::register]]
SEXP agents_sbm_cpp(
    SEXP m,
    const std::vector<int> & block_sizes,
    const std::vector<double> & mixing_matrix,
    bool row_major = true
) {

    // Validating the input (should be all non-negative)
    std::vector<size_t> block_sizes_sz;
    block_sizes_sz.reserve(block_sizes.size());
    for (auto b : block_sizes) {
        if (b < 0)
            stop("Block sizes should be non-negative.\n");

        block_sizes_sz.push_back(static_cast<size_t>(b));
    }

    external_pointer<Model<>> ptr(m);
    ptr->agents_sbm(
        block_sizes_sz,
        mixing_matrix,
        row_major
    );

    return m;

}


[[cpp11::register]]
SEXP agents_from_edgelist_cpp(
  SEXP m,
  const std::vector<int> & source,
  const std::vector<int> & target,
  int size,
  bool directed
) {

  external_pointer<Model<>> ptr(m);
  ptr->agents_from_edgelist(source, target, size, directed);

  return m;

}

[[cpp11::register]]
cpp11::data_frame get_network_cpp(SEXP model) {

  external_pointer<Model<>> modelptr(model);

  std::vector<int> from;
  std::vector<int> to;

  modelptr->write_edgelist(from, to);

  return cpp11::writable::data_frame({
    "from"_nm = from,
    "to"_nm   = to
  });

}
