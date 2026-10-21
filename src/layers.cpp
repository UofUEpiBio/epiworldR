#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

// A layer is held through a std::shared_ptr, so an R layer object and the
// model it is added to share the same layer (a layer created in R) or the R
// object borrows the model's (get_layer(); the R side keeps the model alive).
typedef std::shared_ptr< ContactLayer<> > LayerPtr;
typedef layer::Network<>::Generator Generator;
typedef layer::Network<>::Rewiring Rewiring;

static LayerPtr & layer_ptr(SEXP layer) {
  external_pointer< LayerPtr > ptr(layer);
  return *ptr;
}

static layer::Network<> & network_ref(SEXP layer) {

  auto * net = layer_ptr(layer)->as_network();
  if (net == nullptr)
    stop("The layer is not a network layer.");

  return *net;

}

static layer::Mixing<> & mixing_ref(SEXP layer) {

  auto * mix = dynamic_cast< layer::Mixing<> * >(layer_ptr(layer).get());
  if (mix == nullptr)
    stop("The layer is not a mixing layer.");

  return *mix;

}

// Generators and graphs -------------------------------------------------------

[[cpp11::register]]
SEXP graph_smallworld_cpp(int k, double p, bool directed) {
  return external_pointer< Generator >(
    new Generator(graph::smallworld<>(k, p, directed))
  );
}

[[cpp11::register]]
SEXP graph_bernoulli_cpp(double p, bool directed) {
  return external_pointer< Generator >(
    new Generator(graph::bernoulli<>(p, directed))
  );
}

[[cpp11::register]]
SEXP graph_sbm_cpp(
    const std::vector< int > & block_sizes,
    const std::vector< double > & mixing_matrix,
    bool row_major
) {

  std::vector< size_t > sizes(block_sizes.begin(), block_sizes.end());
  return external_pointer< Generator >(
    new Generator(graph::sbm<>(sizes, mixing_matrix, row_major))
  );

}

[[cpp11::register]]
SEXP graph_edgelist_cpp(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) {
  return external_pointer< AdjList >(
    new AdjList(source, target, size, directed)
  );
}

[[cpp11::register]]
SEXP rewire_degseq_cpp(std::string param) {
  return external_pointer< Rewiring >(
    new Rewiring(graph::rewire_degseq<>(param))
  );
}

// Layers ----------------------------------------------------------------------

[[cpp11::register]]
SEXP layer_network_generator_cpp(std::string name, SEXP generator) {
  external_pointer< Generator > gen(generator);
  return external_pointer< LayerPtr >(
    new LayerPtr(layer::network<>(name, *gen))
  );
}

[[cpp11::register]]
SEXP layer_network_graph_cpp(std::string name, SEXP graph) {
  external_pointer< AdjList > al(graph);
  return external_pointer< LayerPtr >(
    new LayerPtr(layer::network<>(name, *al))
  );
}

[[cpp11::register]]
SEXP add_layer_cpp(SEXP model, SEXP layer) {
  external_pointer< Model<> > ptr(model);
  ptr->add_layer(layer_ptr(layer));
  return model;
}

// Borrows the model's layer: the shared pointer does not own it, so the R
// object must keep the model alive (get_layer() does).
[[cpp11::register]]
SEXP get_layer_cpp(SEXP model, std::string name) {
  external_pointer< Model<> > ptr(model);
  return external_pointer< LayerPtr >(
    new LayerPtr(&ptr->get_layer(name), [](ContactLayer<> *) {})
  );
}

[[cpp11::register]]
std::string get_layer_kind_cpp(SEXP layer) {

  auto * ptr = layer_ptr(layer).get();
  if (ptr->as_network() != nullptr)
    return "network";

  if (dynamic_cast< layer::Mixing<> * >(ptr) != nullptr)
    return "mixing";

  return "other";

}

[[cpp11::register]]
std::string get_layer_name_cpp(SEXP layer) {
  return layer_ptr(layer)->get_name();
}

// Network layers --------------------------------------------------------------

[[cpp11::register]]
bool add_edge_cpp(SEXP layer, int i, int j) {
  return network_ref(layer).add_edge(i, j);
}

[[cpp11::register]]
bool rm_edge_cpp(SEXP layer, int i, int j) {
  return network_ref(layer).rm_edge(i, j);
}

[[cpp11::register]]
bool has_edge_cpp(SEXP layer, int i, int j) {
  return network_ref(layer).has_edge(i, j);
}

[[cpp11::register]]
cpp11::data_frame get_edgelist_cpp(SEXP layer) {

  std::vector< int > from;
  std::vector< int > to;
  network_ref(layer).get_edgelist(from, to);

  return cpp11::writable::data_frame({
    "from"_nm = from,
    "to"_nm   = to
  });

}

[[cpp11::register]]
SEXP set_rewiring_cpp(SEXP layer, SEXP rewiring) {
  external_pointer< Rewiring > fun(rewiring);
  network_ref(layer).set_rewiring(*fun);
  return layer;
}

// Mixing layers ---------------------------------------------------------------

[[cpp11::register]]
std::vector< double > get_contact_matrix_cpp(SEXP layer) {
  return mixing_ref(layer).get_contact_matrix();
}

[[cpp11::register]]
SEXP set_contact_matrix_cpp(SEXP layer, std::vector< double > contact_matrix) {
  mixing_ref(layer).set_contact_matrix(contact_matrix);
  return layer;
}
