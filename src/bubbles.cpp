#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using namespace epiworld;
using namespace cpp11;

[[cpp11::register]]
SEXP bubbles_cpp(
    SEXP model,
    std::vector< int > household_id,
    std::string flavor,
    int group_size,
    double within_factor,
    int start_day,
    int end_day,
    int rewire_every,
    std::string name,
    int max_households
) {

  external_pointer< Model<int> > modelptr(model);

  std::vector< size_t > hh(household_id.begin(), household_id.end());

  BubbleFlavor fl = (flavor == "peer") ?
    BubbleFlavor::Peer : BubbleFlavor::Household;

  Bubbles<int> bubbles(
    hh,
    fl,
    static_cast< size_t >(group_size),
    static_cast< epiworld_double >(within_factor),
    start_day,
    end_day,
    rewire_every,
    name,
    static_cast< size_t >(max_households)
  );

  // deploy() installs the bubble tool and the scheduler event on the model.
  // The intervention's shared state is kept alive by the closures captured in
  // those objects, so the local `bubbles` can safely go out of scope.
  bubbles.deploy(*modelptr);

  return model;

}
