#ifndef EPIWORLDR_ASSIGN_H
#define EPIWORLDR_ASSIGN_H

// An `epiworld_assign` (R/assign.R) as an epiworld::Assign. The R object is a
// list: `form` ("share", "count", "shares", "counts", "group_shares",
// "group_counts", "agents", "labels", "by"), and what the form takes:
// `values` (numbers), `labels` (character or integer), `ids`, `layer`, `fun`.
// A label given as a string is turned into an integer by `label_of` (a stage,
// for viruses and initial stages); an integer label is used as is (a group).

#include <functional>
#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"

using LabelOf = std::function< int(const std::string &) >;

inline int assign_label(SEXP labels, R_xlen_t i, const LabelOf & label_of) {

  if (TYPEOF(labels) == STRSXP)
  {
    if (STRING_ELT(labels, i) == NA_STRING)
      return -1;
    if (!label_of)
      cpp11::stop("This placement takes no labels.");
    return label_of(std::string(CHAR(STRING_ELT(labels, i))));
  }

  // NA (logical) is no label
  if ((TYPEOF(labels) == LGLSXP) && (LOGICAL(labels)[i] == NA_LOGICAL))
    return -1;

  if (TYPEOF(labels) == REALSXP)
    return ISNA(REAL(labels)[i]) ? -1 : static_cast< int >(REAL(labels)[i]);

  if (TYPEOF(labels) != INTSXP)
    cpp11::stop("A label is a string or a whole number.");

  const int l = INTEGER(labels)[i];
  return (l == NA_INTEGER) ? -1 : l;

}

inline epiworld::Assign assign_from_r(SEXP a, const LabelOf & label_of) {

  // No placement: the use's default (e.g., add_virus() without one)
  if (Rf_isNull(a))
    return epiworld::Assign();

  if (!Rf_inherits(a, "epiworld_assign"))
    cpp11::stop("A placement is made by one of the assign_*() functions.");

  cpp11::list al(a);
  const std::string form = cpp11::as_cpp< std::string >(al["form"]);

  if (form == "share")
    return epiworld::assign::share(cpp11::as_cpp< double >(al["values"]));

  if (form == "count")
    return epiworld::assign::count(
      static_cast< size_t >(cpp11::as_cpp< double >(al["values"]))
    );

  if ((form == "shares") || (form == "counts"))
  {
    SEXP labels = al["labels"];
    cpp11::doubles values(al["values"]);
    if (form == "shares")
    {
      std::vector< std::pair< int, epiworld_double > > s;
      for (R_xlen_t i = 0; i < values.size(); ++i)
        s.emplace_back(assign_label(labels, i, label_of), values[i]);
      return epiworld::assign::shares(std::move(s));
    }
    std::vector< std::pair< int, size_t > > c;
    for (R_xlen_t i = 0; i < values.size(); ++i)
      c.emplace_back(
        assign_label(labels, i, label_of),
        static_cast< size_t >(values[i])
      );
    return epiworld::assign::counts(std::move(c));
  }

  if ((form == "group_shares") || (form == "group_counts"))
  {
    const std::string layer = cpp11::as_cpp< std::string >(al["layer"]);
    auto values = cpp11::as_cpp< std::vector< double > >(al["values"]);
    return (form == "group_shares") ?
      epiworld::assign::group_shares(layer, values) :
      epiworld::assign::group_counts(layer, values);
  }

  if (form == "agents")
  {
    std::vector< size_t > ids;
    for (auto i : cpp11::integers(al["ids"]))
      ids.push_back(static_cast< size_t >(i));
    return epiworld::assign::agents(std::move(ids));
  }

  if (form == "labels")
  {
    SEXP labels = al["labels"];
    std::vector< int > v(static_cast< size_t >(Rf_xlength(labels)));
    for (R_xlen_t i = 0; i < Rf_xlength(labels); ++i)
      v[static_cast< size_t >(i)] = assign_label(labels, i, label_of);
    return epiworld::assign::labels(std::move(v));
  }

  if (form == "by")
  {
    // NOTE: calls R, once per agent at each reset; R is not thread-safe, so
    // run_multiple() needs nthreads = 1 (as for globalevent_fun()).
    cpp11::function fun(al["fun"]);
    return epiworld::assign::by(
      [fun, label_of](epiworld::Agent & agent, epiworld::Model & model) -> int {

        cpp11::external_pointer< epiworld::Agent > aptr(&agent, false);
        cpp11::sexp agent_r(aptr);
        agent_r.attr("class") = "epiworld_agent";

        cpp11::external_pointer< epiworld::Model > mptr(&model, false);
        cpp11::sexp model_r(mptr);
        model_r.attr("class") = "epiworld_model";

        cpp11::sexp res = fun(agent_r, model_r);
        if (Rf_length(res) != 1)
          cpp11::stop("assign_by(): the function returns one label per agent.");
        return assign_label(res, 0, label_of);

      }
    );
  }

  cpp11::stop("Unknown placement \"%s\".", form.c_str());

}

#endif
