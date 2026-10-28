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
typedef layer::Mixing<>::Grouping Grouping;

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
SEXP grouping_random_cpp(std::vector< double > sizes, bool as_proportion) {
  return external_pointer< Grouping >(
    new Grouping(grouping::random<>(std::move(sizes), as_proportion))
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
SEXP layer_mixing_cpp(std::string name) {
  return external_pointer< LayerPtr >(new LayerPtr(layer::mixing<>(name)));
}

// Groups from R: 0-based, NA for none (NO_GROUP)
[[cpp11::register]]
SEXP set_groups_cpp(SEXP layer, const std::vector< int > & groups, int n_groups) {

  std::vector< size_t > g(groups.size());
  for (size_t i = 0u; i < groups.size(); ++i)
    g[i] = (groups[i] == NA_INTEGER) ?
      layer::Mixing<>::NO_GROUP : static_cast< size_t >(groups[i]);

  mixing_ref(layer).set_groups(std::move(g), n_groups);
  return layer;

}

[[cpp11::register]]
SEXP set_grouping_cpp(SEXP layer, SEXP grouping) {
  external_pointer< Grouping > gr(grouping);
  mixing_ref(layer).set_groups(*gr);
  return layer;
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
SEXP get_layer_cpp(SEXP model, SEXP which) {

  external_pointer< Model<> > ptr(model);
  auto & layer = (TYPEOF(which) == STRSXP) ?
    ptr->get_layer(as_cpp< std::string >(which)) :
    ptr->get_layer(as_cpp< size_t >(which));

  return external_pointer< LayerPtr >(
    new LayerPtr(&layer, [](ContactLayer<> *) {})
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

[[cpp11::register]]
std::vector< double > get_contact_prob_cpp(SEXP layer) {
  return mixing_ref(layer).get_contact_prob();
}

[[cpp11::register]]
SEXP set_contact_rate_cpp(SEXP layer, std::string param) {
  mixing_ref(layer).set_contact_rate(param);
  return layer;
}

[[cpp11::register]]
SEXP set_contact_prob_cpp(SEXP layer, std::vector< double > contact_prob) {
  mixing_ref(layer).set_contact_prob(contact_prob);
  return layer;
}

[[cpp11::register]]
int get_group_cpp(SEXP layer, int agent) {
  size_t g = mixing_ref(layer).get_group(agent);
  return (g == layer::Mixing<>::NO_GROUP) ? NA_INTEGER : static_cast< int >(g);
}

[[cpp11::register]]
SEXP set_group_cpp(SEXP layer, int agent, int group) {
  mixing_ref(layer).set_group(
    agent,
    (group == NA_INTEGER) ?
      layer::Mixing<>::NO_GROUP : static_cast< size_t >(group)
  );
  return layer;
}

[[cpp11::register]]
std::vector< int > members_of_cpp(SEXP layer, int group) {
  const auto & m = mixing_ref(layer).members_of(group);
  return std::vector< int >(m.begin(), m.end());
}

[[cpp11::register]]
cpp11::data_frame groupmates_of_cpp(SEXP layer, int agent, int from, int to) {

  std::vector< layer::Mixing<>::Groupmate > mates;
  mixing_ref(layer).groupmates_of(agent, from, to, mates);

  std::vector< int > id, last_day, n_days;
  for (const auto & m : mates)
  {
    id.push_back(static_cast< int >(m.agent));
    last_day.push_back(m.last_day);
    n_days.push_back(m.n_days);
  }

  return cpp11::writable::data_frame({
    "agent"_nm    = id,
    "last_day"_nm = last_day,
    "n_days"_nm   = n_days
  });

}

// Every kind of layer -----------------------------------------------------------

[[cpp11::register]]
SEXP set_weight_cpp(SEXP layer, std::string param) {
  layer_ptr(layer)->set_weight(param);
  return layer;
}

[[cpp11::register]]
SEXP add_virus_target_cpp(SEXP layer, int lineage_id) {
  layer_ptr(layer)->add_virus_target(lineage_id);
  return layer;
}

[[cpp11::register]]
SEXP enable_tracing_cpp(SEXP layer, int max_contacts) {
  layer_ptr(layer)->enable_tracing(static_cast< size_t >(max_contacts));
  return layer;
}

[[cpp11::register]]
SEXP disable_tracing_cpp(SEXP layer) {
  layer_ptr(layer)->disable_tracing();
  return layer;
}

[[cpp11::register]]
bool is_tracing_cpp(SEXP layer) {
  return layer_ptr(layer)->is_tracing();
}

[[cpp11::register]]
cpp11::data_frame get_tracing_cpp(SEXP layer) {

  ContactTracing & log = layer_ptr(layer)->get_tracing();

  std::vector< int > source, contact, day;
  for (size_t a = 0u; a < log.get_n_agents(); ++a)
  {
    const size_t n = std::min(log.get_n_contacts(a), log.get_max_contacts());
    for (size_t k = 0u; k < n; ++k)
    {
      auto [c, d] = log.get_contact(a, k);
      source.push_back(static_cast< int >(a));
      contact.push_back(static_cast< int >(c));
      day.push_back(d);
    }
  }

  return cpp11::writable::data_frame({
    "source"_nm  = source,
    "contact"_nm = contact,
    "day"_nm     = day
  });

}
