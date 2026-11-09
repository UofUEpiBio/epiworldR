#ifndef EPIWORLD_MODEL_MEAT_NATURALHISTORY_HPP
#define EPIWORLD_MODEL_MEAT_NATURALHISTORY_HPP

/**
 * @file model-meat-naturalhistory.hpp
 * @brief The model's side of natural histories and statuses.
 *
 * @details An agent has a stage (from its virus's natural history) and a
 * status (from the model's `Statuses`). The engine keeps one index for the
 * pair, a *cell*: `cell = stage * n_statuses + status`, which is the agent's
 * `state`. With one status, the cells are the stages. Events carry the target
 * stage or status; the cell is computed when the event is applied, so a stage
 * change and a status change on the same day combine.
 */

// Installing ------------------------------------------------------------------

inline void Model::nh_install(std::shared_ptr< NaturalHistory > def)
{

    if (nstates != 0u)
        throw std::logic_error(
            "The model already has states. A virus with a natural history "
            "gives the model its stages as its states: add it before anything "
            "else."
        );

    if (def->is_carrier(0))
        throw std::logic_error(
            "The first stage of \"" + def->get_name() + "\" (\"" +
            def->get_stages()[0] + "\") is where agents start, so it cannot be "
            "a carrier stage."
        );

    bool default_infection = false;
    bool infectious = false;
    for (int s = 0; s < static_cast< int >(def->get_n_stages()); ++s)
    {
        const auto & info = def->get_stage_info(s);
        if ((info.infection >= 0) && !def->get_transition(info.infection).sampler)
            default_infection = true;
        if (def->is_infectious(s))
            infectious = true;
    }

    if (default_infection && !infectious)
        throw std::logic_error(
            "The natural history \"" + def->get_name() + "\" has an infection "
            "by contact but no infectious stage (stage::infectious)."
        );

    nh = std::move(def);
    nh_build_states();

}

inline void Model::set_statuses(Statuses s)
{

    if (running)
        throw std::logic_error("set_statuses() is called between runs.");

    if (s.get_n_stages() == 0u)
        throw std::logic_error("The statuses have no status.");

    s.finalize(*this, "Statuses");
    st = std::make_shared< const Statuses >(std::move(s));

    if (nh)
        nh_build_states();

}

inline void Model::nh_build_states()
{

    const int ns = static_cast< int >(nh->get_n_stages());
    const int nk = static_cast< int >(nh_nstatus());

    // A status transition's stages (`.when()` and a duration's `.since()`)
    // are the virus's, and its name is not one of the virus's
    // (transition_id() finds either by name)
    if (st)
        for (size_t t = 0u; t < st->get_n_transitions(); ++t)
        {
            const auto & tr = st->get_transition(static_cast< int >(t));
            for (int s : tr.when)
                nh->check_stage(s, "when");
            if (
                (tr.kind == TransitionKind::after) &&
                (tr.duration.get_anchor() >= 0)
            )
                nh->check_stage(tr.duration.get_anchor(), "anchor");
            for (size_t u = 0u; u < nh->get_n_transitions(); ++u)
                if (nh->get_transition(static_cast< int >(u)).name == tr.name)
                    throw std::logic_error(
                        "The virus and the statuses both have a transition "
                        "named \"" + tr.name + "\"; give them different names."
                    );
        }

    // Every cell with a drawn exit (of its stage or its status) is updated
    UpdateFun update = [](Agent * p, Model * m) -> void {
        m->nh_step(*p);
    };

    states_labels.clear();
    state_fun.clear();
    nstates = 0u;
    infectious_states.clear();
    excluded_states.clear();
    nh_queue.assign(static_cast< size_t >(ns * nk), 0);

    for (int s = 0; s < ns; ++s)
        for (int k = 0; k < nk; ++k)
        {

            const auto & sinfo = nh->get_stage_info(s);
            const bool status_drawn = st && st->get_stage_info(k).needs_visit;
            const bool drawn = (sinfo.infection >= 0) || sinfo.needs_visit ||
                status_drawn;

            states_labels.push_back(
                (nk == 1) ? nh->get_stages()[s] :
                st->get_stages()[k] + " / " + nh->get_stages()[s]
            );
            state_fun.push_back(drawn ? update : nullptr);

            const int c = s * nk + k;
            const bool excluded = (nh->get_contact(s) <= 0.0) ||
                (st && (st->get_contact(k) <= 0.0));

            if (excluded)
                excluded_states.push_back(c);
            else if (nh->is_infectious(s))
                infectious_states.push_back(c);

            // Carriers in the layers queue their contacts; the others with a
            // drawn exit queue themselves, except those whose only exit is an
            // infection drawn from their contacts (the default sampler),
            // which the carriers queue. A custom sampler need not use the
            // contacts, so its agents queue themselves.
            const bool contact_driven = sinfo.pushable && !status_drawn;
            nh_queue[c] = (nh->is_carrier(s) && !excluded) ? 2 :
                ((drawn && !contact_driven) ? 1 : 0);

        }

    nstates = static_cast< epiworld_fast_uint >(ns * nk);

    // The index has one slot per state; it is rebuilt when the next run starts
    state_index_ready = false;

}

// Accessors -------------------------------------------------------------------

inline int Model::get_stage(const Agent & agent) const
{
    if (!nh)
        throw std::logic_error("The model has no natural history.");
    return nh_stage_of(agent.state);
}

inline int Model::get_status(const Agent & agent) const
{
    if (!nh)
        throw std::logic_error("The model has no natural history.");
    return nh_status_of(agent.state);
}

inline int Model::transition_id(std::string_view name) const
{

    if (!nh)
        throw std::logic_error("The model has no natural history.");

    for (size_t t = 0u; t < nh->get_n_transitions(); ++t)
        if (nh->get_transition(static_cast< int >(t)).name == name)
            return static_cast< int >(t);

    if (st)
        return st->transition_id(name);

    throw std::logic_error(
        "The model has no transition \"" + std::string(name) + "\"."
    );

}

inline int Model::get_stage_entered(size_t id, int stage) const
{
    if (!nh)
        throw std::logic_error("The model has no natural history.");
    nh->check_stage(stage, "");
    const int row = nh_row.at(id);
    if (row < 0)
        return (stage == 0) ? 0 : -1;
    return nh_rows[static_cast< size_t >(row) * (1u + nh->get_n_stages()) +
        1u + static_cast< size_t >(stage)];
}

inline int Model::get_stage_until(size_t id) const
{
    if (!nh)
        throw std::logic_error("The model has no natural history.");
    return nh_until.at(id);
}

inline int Model::get_episode(size_t id) const
{
    if (!nh)
        throw std::logic_error("The model has no natural history.");
    const int row = nh_row.at(id);
    if (row < 0)
        return 0;
    return nh_rows[static_cast< size_t >(row) * (1u + nh->get_n_stages())];
}

inline int * Model::nh_record(size_t id)
{

    const size_t width = 1u + nh->get_n_stages();
    if (nh_row[id] < 0)
    {
        // No infection yet, in the first stage since day 0
        nh_row[id] = static_cast< int >(nh_rows.size() / width);
        nh_rows.resize(nh_rows.size() + width, -1);
        nh_rows[nh_rows.size() - width]      = 0;
        nh_rows[nh_rows.size() - width + 1u] = 0;
    }

    return nh_rows.data() + static_cast< size_t >(nh_row[id]) * width;

}

inline bool Model::nh_removed(int stage) const
{
    // A stage that takes agents out of every layer and has no exit
    if (nh->get_contact(stage) > 0.0)
        return false;
    for (size_t t = 0u; t < nh->get_n_transitions(); ++t)
        if (nh->is_source(static_cast< int >(t), stage))
            return false;
    return true;
}

// Moving agents from code ------------------------------------------------------

inline bool Model::fire(int transition, Agent & agent)
{

    if (!nh)
        throw std::logic_error("fire() needs a natural history.");

    const int s = nh_stage_of(agent.state);

    if (nh_is_status_tr(transition))
    {

        const int k = nh_status_tr(transition);
        if (!st || (static_cast< size_t >(k) >= st->get_n_transitions()))
            throw std::range_error(
                "There is no status transition " + std::to_string(k) + "."
            );

        if (!st->is_source(k, nh_status_of(agent.state)) || nh_removed(s))
            return false;

    }
    else
    {

        if ((transition < 0) || (static_cast< size_t >(transition) >= nh->get_n_transitions()))
            throw std::range_error(
                "There is no transition " + std::to_string(transition) + "."
            );

        if (nh->get_transition(transition).kind == TransitionKind::infection)
            throw std::logic_error(
                "An infection is not fired; seed the virus instead."
            );

        if (!nh->is_source(transition, s))
            return false;

    }

    nh_propose(transition, agent);
    return true;

}

inline void Model::react(Trigger on, int transition)
{
    if (!nh)
        throw std::logic_error("react() needs a natural history.");
    nh_reactions.push_back({std::move(on), transition, nullptr});
}

inline void Model::react(
    Trigger on,
    std::function< void(Agent &, Model &) > fun
)
{
    if (!nh)
        throw std::logic_error("react() needs a natural history.");
    nh_reactions.push_back({std::move(on), -1, std::move(fun)});
}

inline void Model::set_initial(Assign a)
{

    if (!nh)
        throw std::logic_error("set_initial() needs a natural history.");

    nh_initial = std::move(a);

}

// The run ---------------------------------------------------------------------

inline void Model::nh_reset()
{

    const size_t n = population.size();

    nh_until.assign(n, INT_MAX);
    nh_row.assign(n, -1);
    nh_rows.clear();
    nh_status_until.assign(st ? n : 0u, INT_MAX);
    nh_status_sod.assign(st ? n : 0u, 0);
    nh_pending = -1;

    // Everyone starts in the first stage and the first status, on day 0
    const int after = nh->stages[0u].after;
    const int s_after = st ? st->stages[0u].after : -1;
    if ((after >= 0) || (s_after >= 0))
        for (auto & p : population)
        {
            const size_t id = static_cast< size_t >(p.id);
            if (after >= 0)
                nh_until[id] = nh->transitions[after].duration.draw(p, *this);
            if (s_after >= 0)
                nh_status_until[id] = st->transitions[s_after].duration.draw(p, *this);
        }

    if (use_queuing && (nh_queue[0u] == 1))
        for (size_t id = 0u; id < n; ++id)
            queue.shift(id, 1);

}

inline void Model::nh_step(Agent & p)
{

    const unsigned int c = p.state;
    const int s = nh_stage_of(c);
    const StageInfo & si = nh->stages[static_cast< size_t >(s)];

    // The stage's proposal --------------------------------------------------
    bool proposed = false;

    // The values that apply to the agent (its virus's, or the model's)
    const NaturalHistory & h = nh_of(p);

    // 1. Infection by contact (not out of the layers)
    if (
        (si.infection >= 0) && (p.virus == nullptr) &&
        (state_excluded.empty() || !state_excluded[c])
    )
    {

        const auto & tr = h.transitions[si.infection];
        Virus * v = tr.sampler ? tr.sampler(p, *this) : sample_virus(&p);
        if (v != nullptr)
        {
            nh_propose(si.infection, p, v);
            proposed = true;
        }

    }

    // 2. The rates compete: one draw picks one, or none
    const size_t nr = si.rates.size();
    if (!proposed && (nr > 0u))
    {

        reserve_tmp(nr);
        for (size_t k = 0u; k < nr; ++k)
        {

            const auto & tr = h.transitions[si.rates[k]];
            epiworld_double x = tr.value.eval(p, *this);

            if ((tr.effect != Effect::none) && (p.virus != nullptr))
            {
                if (tr.effect == Effect::recovery)
                    x = 1.0 - (1.0 - x) * (1.0 - p.get_recovery_enhancer(p.virus, *this));
                else
                    x = x * (1.0 - p.get_death_reduction(p.virus, *this));
            }

            array_double_tmp[k] = x;

        }

        const int which = roulette(static_cast< epiworld_fast_uint >(nr), this);
        if (which >= 0)
        {
            nh_propose(si.rates[static_cast< size_t >(which)], p);
            proposed = true;
        }

    }

    // 3. The duration, if it ended (a deadline: the rates went first)
    if (!proposed && (si.after >= 0) && (nh_until[static_cast< size_t >(p.id)] <= today()))
        nh_propose(si.after, p);

    // The status's proposal (judged on the stage at the start of the day)
    if (st)
        nh_step_status(p, c, s);

}

inline void Model::nh_step_status(Agent & p, unsigned int c, int s)
{

    const int k = nh_status_of(c);
    const StageInfo & ki = st->stages[static_cast< size_t >(k)];

    auto applies = [s](const Transition & tr) -> bool {
        return tr.when.empty() ||
            (std::find(tr.when.begin(), tr.when.end(), s) != tr.when.end());
    };

    const size_t nkr = ki.rates.size();
    if (nkr > 0u)
    {

        reserve_tmp(nkr);
        for (size_t r = 0u; r < nkr; ++r)
        {
            const auto & tr = st->transitions[ki.rates[r]];
            array_double_tmp[r] = applies(tr) ? tr.value.eval(p, *this) : 0.0;
        }

        const int which = roulette(static_cast< epiworld_fast_uint >(nkr), this);
        if (which >= 0)
        {
            nh_propose(-2 - ki.rates[static_cast< size_t >(which)], p);
            return;
        }

    }

    if (
        (ki.after >= 0) &&
        (nh_status_until[static_cast< size_t >(p.id)] <= today()) &&
        applies(st->transitions[ki.after])
    )
        nh_propose(-2 - ki.after, p);

}

inline void Model::nh_propose(
    int transition,
    Agent & p,
    Virus * virus
)
{

    const bool status = nh_is_status_tr(transition);

    // The infection's target is the infecting virus's; the others, the agent's
    // (a status transition, the model's statuses)
    const auto & tr = status ?
        st->transitions[nh_status_tr(transition)] :
        (((virus != nullptr) && virus->nh_def) ?
            virus->nh_def->transitions[transition] :
            nh_of(p).transitions[transition]);

    // A branch draws its target
    int to = tr.to;
    if (!tr.branch_to.empty())
    {

        const epiworld_double u = runif();
        epiworld_double cum = 0.0;
        for (size_t k = 0u; k < tr.branch_to.size(); ++k)
        {
            cum += tr.branch_prob[k].eval(p, *this);
            to = tr.branch_to[k];
            if (u < cum)
                break;
        }

    }

    // The event carries the target (a stage or a status); events_run() makes
    // it a cell when it applies it
    nh_pending = transition;

    const int s = nh_stage_of(p.state);
    if (status)
        p.change_state(*this, static_cast< epiworld_fast_uint >(to), Queue::NoOne);
    else if (tr.kind == TransitionKind::infection)
        p.set_virus(*this, *virus, to, Queue::NoOne);
    else if (nh->is_carrier(s) && !nh->is_carrier(to))
        p.rm_virus(*this, to, Queue::NoOne);
    else
        p.change_state(*this, static_cast< epiworld_fast_uint >(to), Queue::NoOne);

    nh_pending = -1;

}

inline void Model::nh_enter(
    Agent & p,
    unsigned int stage_old,
    int transition
)
{

    const size_t id = static_cast< size_t >(p.id);
    const size_t ns = nh->get_n_stages();
    const int s = nh_stage_of(p.state);

    // A new infection starts a new episode: its entry days start afresh
    const bool infected =
        (transition == NH_SEED) ||
        (!nh->is_carrier(static_cast< int >(stage_old)) && nh->is_carrier(s));

    int * record  = nh_record(id);
    int * entered = record + 1;
    if (infected)
    {
        ++record[0];
        std::fill(entered, entered + ns, -1);
    }

    entered[s] = today();

    const int after = nh->stages[static_cast< size_t >(s)].after;
    if (after >= 0)
    {

        const auto & d = nh_of(p).transitions[after].duration;
        int base = today();
        if ((d.get_anchor() >= 0) && (entered[d.get_anchor()] >= 0))
            base = entered[d.get_anchor()];

        nh_until[id] = base + d.draw(p, *this);

    }
    else
        nh_until[id] = INT_MAX;

}

inline void Model::nh_react(Agent & p, int transition, int stage_old)
{

    if (nh_reactions.empty())
        return;

    const size_t id = static_cast< size_t >(p.id);
    const int s = nh_stage_of(p.state);

    // The status at the start of the day (the agent changed today)
    const int status_sod = st ? nh_status_sod[id] : 0;

    for (const auto & r : nh_reactions)
    {

        const bool hit = (r.on.kind == Trigger::Kind::entered) ?
            ((s == r.on.target) && (stage_old != s)) :
            (transition == r.on.target);

        if (!hit)
            continue;

        if (
            !r.on.was_status.empty() &&
            (std::find(r.on.was_status.begin(), r.on.was_status.end(), status_sod) ==
                r.on.was_status.end())
        )
            continue;

        if (r.fun)
            r.fun(p, *this);
        else
            fire(r.transition, p);

    }

}

inline void Model::nh_queue_move(Agent & p, unsigned int state_old)
{

    const char from = nh_queue[state_old];
    const char to = nh_queue[p.state];
    if (from == to)
        return;

    const size_t id = static_cast< size_t >(p.id);

    if (from == 2)
        queue -= &p;
    else if (from == 1)
        queue.shift(id, -1);

    if (to == 2)
        queue += &p;
    else if (to == 1)
        queue.shift(id, 1);

}

inline void Model::nh_initial_place()
{

    if (!nh_initial.is_set())
        return;

    // The agents left without the virus, in the first stage
    auto & buf = assign_buffers;
    buf.eligible.resize(population.size());
    for (auto & p : population)
        buf.eligible[static_cast< size_t >(p.id)] =
            (p.virus == nullptr) && (nh_stage_of(p.state) == 0);

    nh_initial.draw(*this, buf);

    for (const auto & x : buf.out)
    {

        const size_t i = x.first;
        const int s    = x.second;

        if (
            (s < 0) || (static_cast< size_t >(s) >= nh->get_n_stages()) ||
            nh->is_carrier(s)
        )
            throw std::logic_error(
                "set_initial(): the label " + std::to_string(s) + " is not a "
                "stage without the virus (an initial stage is a stage, e.g., "
                "assign::shares({{R, 0.1}}); carriers are seeded with the "
                "virus)."
            );

        nh_pending = NH_INITIAL;
        population[i].change_state(
            *this, static_cast< epiworld_fast_uint >(s), Queue::NoOne
        );
        nh_pending = -1;

    }

    events_run();

}

inline epiworld_double Model::nh_infection_prob(
    Agent & target,
    const Virus & v
)
{
    const int tr = nh->stages[static_cast< size_t >(nh_stage_of(target.state))].infection;
    if (tr < 0)
        return 0.0;

    // The infector's virus says how infectious it is
    const auto & h = v.nh_def ? *v.nh_def : *nh;
    return h.transitions[tr].value.eval(target, &v, *this);
}

inline std::string Model::draw_natural_history() const
{

    if (!nh)
        throw std::logic_error(
            "The model has no natural history to draw (add a virus first)."
        );

    std::string out = "flowchart LR\n";

    auto subgraph = [&out](
        const std::string & id,
        const std::string & title,
        const std::string & body
    ) {
        out += "    subgraph " + id + "[\"" + title + "\"]\n";
        out += "        direction LR\n";
        std::istringstream lines(body);
        std::string line;
        while (std::getline(lines, line))
            out += "    " + line + "\n";
        out += "    end\n";
    };

    subgraph("stages", "Stages: " + nh->get_name(), nh->draw_mermaid("s"));
    if (st)
        subgraph("statuses", "Statuses", st->draw_mermaid("q", nh.get()));

    out +=
        "    classDef infectious fill:#f8d7da,stroke:#b02a37\n"
        "    classDef carrier fill:#fff3cd,stroke:#997404\n"
        "    classDef away fill:#e2e3e5,stroke:#6c757d\n";

    return out;

}

#endif
