#ifndef EPIWORLDR_NATURAL_HISTORY_H
#define EPIWORLDR_NATURAL_HISTORY_H

// Conversions shared by the natural-history bindings (virus.cpp): the values
// and durations of the transitions, and stages named by their labels.

#include <functional>
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

// A transition's value from R: a number, a model parameter (its name), or an
// `epiworld_value` (value_inverse(), value_complement(), value_fun(), ...)
inline epiworld::Value value_from_r(SEXP x) {

  if (Rf_inherits(x, "epiworld_value"))
    return *cpp11::external_pointer<epiworld::Value>(x);

  if ((TYPEOF(x) == STRSXP) && (Rf_length(x) == 1))
    return epiworld::Value(std::string(CHAR(STRING_ELT(x, 0))));

  if (((TYPEOF(x) == REALSXP) || (TYPEOF(x) == INTSXP)) && (Rf_length(x) == 1))
    return epiworld::Value(Rf_asReal(x));

  cpp11::stop(
    "A value is a number, the name of a model parameter, or an object made "
    "by value_inverse(), value_complement() or value_fun()."
  );

}

// The natural history of `virus` (throws if it has no stages yet)
inline const epiworld::NaturalHistory & virus_nh(const epiworld::Virus & virus) {

  const epiworld::NaturalHistory * nh = virus.get_natural_history();
  if (nh == nullptr)
    cpp11::stop(
      "The virus \"%s\" has no stages yet: add them with add_stage().",
      virus.get_name().c_str()
    );

  return *nh;

}

// The stage of `virus` labeled `label` (throws if none)
inline int stage_from_label(const epiworld::Virus & virus, const std::string & label) {
  return virus_nh(virus).stage_id(label);
}

inline epiworld::Nodes stages_from_labels(
  const epiworld::Virus & virus,
  const cpp11::strings & labels
) {

  std::vector< int > ids;
  for (const auto & l : labels)
    ids.push_back(stage_from_label(virus, std::string(l)));

  return epiworld::Nodes(ids);

}

// A placement's stage labels, resolved against a copy of the stages' labels.
// The resolver is kept with the placement (assign_by() calls it at every
// reset, and a cloned model copies it), so it owns what it reads: it never
// refers to the caller's virus or to the model it was made from.
inline std::function< int(const std::string &) > stage_resolver(
  const epiworld::NaturalHistory & nh
) {

  return [name = nh.get_name(), labels = nh.get_stages()](
    const std::string & label
  ) -> int {
    for (size_t s = 0u; s < labels.size(); ++s)
      if (labels[s] == label)
        return static_cast< int >(s);
    cpp11::stop(
      "The natural history \"%s\" has no stage \"%s\".",
      name.c_str(), label.c_str()
    );
  };

}

#endif
