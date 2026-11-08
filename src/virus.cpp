#include "cpp11.hpp"
#include "cpp11/external_pointer.hpp"
#include "epiworld-common.h"
#include "natural-history.h"
#include "assign.h"

using namespace cpp11;
using namespace epiworld;

#define WrapVirus(a) external_pointer<Virus> (a)

[[cpp11::register]]
SEXP virus_cpp(std::string name) {

  WrapVirus(virus)(new Virus(name));
  return virus;

}

// The natural history ---------------------------------------------------------

[[cpp11::register]]
int add_stage_cpp(SEXP virus, std::string label, bool carrier, bool infectious) {

  unsigned int flags = 0u;
  if (infectious)
    flags = stage::infectious;
  else if (carrier)
    flags = stage::carrier;

  return external_pointer<Virus>(virus)->add_stage(label, flags);

}

[[cpp11::register]]
SEXP set_contact_cpp(SEXP virus, std::string stage, double weight) {

  WrapVirus(vptr)(virus);
  vptr->set_contact(stage_from_label(*vptr, stage), weight);
  return virus;

}

// Names the transition just declared, if given a name
static int name_transition(TransitionRef ref, SEXP name) {

  if (!Rf_isNull(name))
    ref.name(cpp11::as_cpp< std::string >(name));

  return ref.get_id();

}

[[cpp11::register]]
int add_infection_cpp(
  SEXP virus, strings from, std::string to, SEXP rate, SEXP name
) {

  WrapVirus(vptr)(virus);
  return name_transition(
    vptr->infection(
      stages_from_labels(*vptr, from),
      stage_from_label(*vptr, to),
      value_from_r(rate)
    ),
    name
  );

}

[[cpp11::register]]
int add_rate_cpp(
  SEXP virus, strings from, std::string to, SEXP rate, std::string effect,
  SEXP name
) {

  WrapVirus(vptr)(virus);
  TransitionRef ref = vptr->rate(
    stages_from_labels(*vptr, from),
    stage_from_label(*vptr, to),
    value_from_r(rate)
  );

  if (effect == "recovery")
    ref.effect(Effect::recovery);
  else if (effect == "death")
    ref.effect(Effect::death);
  else if (effect != "none")
    cpp11::stop("The effect is one of \"none\", \"recovery\" and \"death\".");

  return name_transition(ref, name);

}

[[cpp11::register]]
int add_after_cpp(
  SEXP virus, std::string from, std::string to, SEXP duration, SEXP since,
  SEXP name
) {

  WrapVirus(vptr)(virus);
  Duration d = *external_pointer<Duration>(duration);
  if (!Rf_isNull(since))
    d = d.since(stage_from_label(*vptr, cpp11::as_cpp< std::string >(since)));

  return name_transition(
    vptr->after(stage_from_label(*vptr, from), stage_from_label(*vptr, to), d),
    name
  );

}

[[cpp11::register]]
int add_external_cpp(SEXP virus, strings from, std::string to, std::string name) {

  WrapVirus(vptr)(virus);
  return vptr->external(
    stages_from_labels(*vptr, from), stage_from_label(*vptr, to), name
  ).get_id();

}

[[cpp11::register]]
SEXP use_param_cpp(SEXP virus, std::string param, std::string replacement) {

  external_pointer<Virus>(virus)->use_param(param, replacement);
  return virus;

}

[[cpp11::register]]
cpp11::writable::list get_stages_virus_cpp(SEXP virus) {

  WrapVirus(vptr)(virus);
  const NaturalHistory * nh = vptr->get_natural_history();

  cpp11::writable::strings labels;
  cpp11::writable::logicals carrier, infectious;
  cpp11::writable::doubles contact;
  if (nh != nullptr)
    for (int s = 0; s < static_cast< int >(nh->get_n_stages()); ++s)
    {
      labels.push_back(nh->get_stages()[static_cast< size_t >(s)]);
      carrier.push_back(nh->is_carrier(s) ? TRUE : FALSE);
      infectious.push_back(nh->is_infectious(s) ? TRUE : FALSE);
      contact.push_back(nh->get_contact(s));
    }

  return cpp11::writable::list({
    "stage"_nm = labels,
    "carrier"_nm = carrier,
    "infectious"_nm = infectious,
    "contact"_nm = contact
  });

}

// Values and durations --------------------------------------------------------

static SEXP wrap_value(Value v) {

  external_pointer<Value> ptr(new Value(std::move(v)));
  cpp11::sexp res(ptr);
  res.attr("class") = "epiworld_value";
  return res;

}

[[cpp11::register]]
SEXP value_inverse_cpp(std::string param) {
  return wrap_value(Value::inverse(param));
}

[[cpp11::register]]
SEXP value_complement_cpp(std::string param) {
  return wrap_value(Value::complement(param));
}

[[cpp11::register]]
SEXP value_fun_cpp(cpp11::function fun) {

  // NOTE: calls R at every evaluation (once per agent and day it applies to).
  // R is not thread-safe: run_multiple() needs nthreads = 1 (as for
  // globalevent_fun()).
  return wrap_value(Value(ValueFun(
    [fun](Agent & agent, Model & model) -> epiworld_double {

      external_pointer<Agent> aptr(&agent, false);
      cpp11::sexp agent_r(aptr);
      agent_r.attr("class") = "epiworld_agent";

      external_pointer<Model> mptr(&model, false);
      cpp11::sexp model_r(mptr);
      model_r.attr("class") = "epiworld_model";

      return cpp11::as_cpp< double >(fun(agent_r, model_r));

    }
  )));

}

[[cpp11::register]]
SEXP rate_fun_logit_cpp(integers vars, doubles coefs, SEXP model) {

  std::vector< int > v = as_cpp< std::vector< int > >(vars);
  std::vector< double > b = as_cpp< std::vector< double > >(coefs);

  if (b.empty())
    cpp11::stop("rate_fun_logit(): -coefs- has at least one element.");

  if (b.size() != v.size())
    cpp11::stop("rate_fun_logit(): -vars- and -coefs- have the same length.");

  const int K = static_cast< int >(
    external_pointer<Model>(model)->get_agents_data_ncols()
  );
  for (auto k : v)
    if ((k < 0) || (k >= K))
      cpp11::stop(
        "rate_fun_logit(): the variable %d is out of range (the agents have "
        "%d).", k, K
      );

  return wrap_value(Value(ValueFun(
    [v, b](Agent & agent, Model & model) -> epiworld_double {

      epiworld_double eta = 0.0;
      for (size_t k = 0u; k < v.size(); ++k)
        eta += agent(static_cast< size_t >(v[k]), model) * b[k];

      return 1.0 / (1.0 + std::exp(-eta));

    }
  )));

}

[[cpp11::register]]
SEXP duration_cpp(std::string kind, SEXP mean, SEXP sd) {

  Duration d;
  if (kind == "fixed")
    d = Duration::fixed(value_from_r(mean));
  else if (kind == "geometric")
    d = Duration::geometric(value_from_r(mean));
  else if (kind == "lognormal")
    d = Duration::lognormal(value_from_r(mean), value_from_r(sd));
  else if (kind == "gamma")
    d = Duration::gamma(value_from_r(mean), value_from_r(sd));
  else
    cpp11::stop("Unknown duration \"%s\".", kind.c_str());

  external_pointer<Duration> ptr(new Duration(d));
  cpp11::sexp res(ptr);
  res.attr("class") = "epiworld_duration";
  return res;

}

// In the model ----------------------------------------------------------------

[[cpp11::register]]
SEXP add_virus_cpp(SEXP m, SEXP v, SEXP assign) {

  WrapVirus(vptr)(v);
  const Virus & virus = *vptr;

  // A placement's labels are the virus's stages
  external_pointer<Model>(m)->add_virus(
    virus,
    assign_from_r(assign, [&virus](const std::string & label) {
      return stage_from_label(virus, label);
    })
  );

  return m;

}

[[cpp11::register]]
SEXP assign_virus_cpp(SEXP m, int virus_pos, SEXP assign) {

  external_pointer<Model> mptr(m);
  const Virus & virus = mptr->get_virus(static_cast< size_t >(virus_pos));
  mptr->assign_virus(
    static_cast< size_t >(virus_pos),
    assign_from_r(assign, [&virus](const std::string & label) {
      return stage_from_label(virus, label);
    })
  );

  return m;

}

[[cpp11::register]]
SEXP rm_virus_cpp(SEXP m, size_t virus_pos) {

  external_pointer<Model>(m)->rm_virus(virus_pos);
  return m;

}

[[cpp11::register]]
SEXP print_virus_cpp(SEXP v) {

  WrapVirus(vptr)(v);
  vptr->print();
  return v;

}

[[cpp11::register]]
std::string get_name_virus_cpp(SEXP virus) {
  return external_pointer<Virus>(virus)->get_name();
}

[[cpp11::register]]
int get_lineage_virus_cpp(SEXP virus) {
  return external_pointer<Virus>(virus)->get_lineage_id();
}

[[cpp11::register]]
SEXP set_name_virus_cpp(SEXP virus, std::string name) {
  external_pointer<Virus>(virus)->set_name(name);
  return virus;
}

[[cpp11::register]]
SEXP clone_virus_cpp(SEXP virus) {

  // clone_ptr() keeps the virus's type. The copy is a new variant: the model
  // gives it its own id when it is added.
  auto cloned = external_pointer<Virus>(virus)->clone_ptr();
  cloned->set_id(-99);

  external_pointer<Virus> ptr(cloned.release());
  return ptr;

}

#undef WrapVirus
