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

    if (tracks.empty() && (nstates != 0u))
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

    tracks.emplace_back();
    tracks.back().id = static_cast< int >(tracks.size() - 1u);
    tracks.back().def = std::move(def);
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

    if (!tracks.empty())
        nh_build_states();

}

/// The most cells (status × each disease's stage) a model can have
constexpr size_t EPI_MAX_CELLS = size_t(1u) << 16;

inline void Model::nh_build_states()
{

    const size_t nk = nh_nstatus();
    const bool several = tracks.size() > 1u;

    // Whether a disease removes agents (a stage out of every layer, without
    // exits); then the others have "Removed (other disease)"
    auto removes = [](const NaturalHistory & h) -> bool {
        for (int s = 0; s < static_cast< int >(h.get_n_stages()); ++s)
        {
            if (h.get_contact(s) > 0.0)
                continue;
            bool exit = false;
            for (size_t t = 0u; t < h.get_n_transitions(); ++t)
                exit = exit || h.is_source(static_cast< int >(t), s);
            if (!exit)
                return true;
        }
        return false;
    };

    size_t ncells = nk;
    for (auto & track : tracks)
    {

        bool other_removes = false;
        for (const auto & other : tracks)
            if ((other.id != track.id) && removes(*other.def))
                other_removes = true;

        if (several && other_removes)
        {
            auto nh = std::make_shared< NaturalHistory >(*track.def);
            nh->stage_labels.push_back("Removed (other disease)");
            nh->stage_flags.push_back(0u);
            nh->stage_contact.push_back(0.0);
            nh->stages.push_back(StageInfo());
            track.removed_other = static_cast< int >(nh->get_n_stages()) - 1;
            track.nh = std::move(nh);
        }
        else
        {
            track.nh = track.def;
            track.removed_other = -1;
        }

        track.n = static_cast< unsigned int >(track.nh->get_n_stages());
        track.stride = static_cast< unsigned int >(ncells);
        ncells *= track.n;

        if (ncells > EPI_MAX_CELLS)
            throw std::logic_error(
                "The model would have " + std::to_string(ncells) + "+ states "
                "(a status and a stage in each disease, combined); the most "
                "is " + std::to_string(EPI_MAX_CELLS) + ". Use fewer diseases, "
                "stages or statuses."
            );

    }

    // A status transition's name is not one of the diseases' (transition_id()
    // finds either by name)
    if (st)
        for (size_t t = 0u; t < st->get_n_transitions(); ++t)
        {
            const auto & tr = st->get_transition(static_cast< int >(t));
            for (const auto & track : tracks)
                for (size_t u = 0u; u < track.nh->get_n_transitions(); ++u)
                    if (track.nh->get_transition(static_cast< int >(u)).name == tr.name)
                        throw std::logic_error(
                            "The virus and the statuses both have a transition "
                            "named \"" + tr.name + "\"; give them different names."
                        );
        }

    // Every cell with a drawn exit (of a stage or of the status) is updated
    // (one disease: its stage, infection first, then the status)
    UpdateFun update = several ?
        UpdateFun([](Agent * p, Model * m) -> void { m->nh_step(*p); }) :
        UpdateFun([](Agent * p, Model * m) -> void {
            const unsigned int c = p->state;
            m->nh_step_stage(m->tracks[0u], *p, c, true);
            if (m->st)
                m->nh_step_status(*p, c);
        });

    for (auto & track : tracks)
    {
        track.stage_at.resize(ncells);
        for (size_t c = 0u; c < ncells; ++c)
            track.stage_at[c] = static_cast< int >((c / track.stride) % track.n);
    }

    states_labels.clear();
    state_fun.clear();
    infectious_states.clear();
    excluded_states.clear();
    nh_queue.assign(ncells, 0);

    for (size_t c = 0u; c < ncells; ++c)
    {

        const int k = static_cast< int >(c % nk);
        const bool status_drawn = st && st->get_stage_info(k).needs_visit;

        bool drawn = status_drawn;
        bool self = status_drawn;   // A reason to visit that is not a contact
        bool excluded = st && (st->get_contact(k) <= 0.0);
        bool infectious = false;
        bool carrier = false;
        std::string label = (nk == 1u) ? std::string() : st->get_stages()[static_cast< size_t >(k)];

        for (const auto & track : tracks)
        {

            const int s = static_cast< int >((c / track.stride) % track.n);
            const auto & si = track.nh->get_stage_info(s);
            const bool stage_drawn = (si.infection >= 0) || si.needs_visit;

            drawn = drawn || stage_drawn;
            // A stage whose only drawn exit is a default infection is visited
            // when a carrier queues it
            self = self || (stage_drawn && !si.pushable);
            excluded = excluded || (track.nh->get_contact(s) <= 0.0);
            infectious = infectious || track.nh->is_infectious(s);
            carrier = carrier || track.nh->is_carrier(s);
            label += (label.empty() ? "" : " / ") + track.nh->get_stages()[static_cast< size_t >(s)];

        }

        states_labels.push_back(std::move(label));
        state_fun.push_back(drawn ? update : nullptr);

        if (excluded)
            excluded_states.push_back(static_cast< int >(c));
        else if (infectious)
            infectious_states.push_back(static_cast< int >(c));

        // Carriers in the layers queue their contacts; the others with a
        // reason to be visited queue themselves
        nh_queue[c] = (carrier && !excluded) ? 2 : (self ? 1 : 0);

    }

    nstates = static_cast< epiworld_fast_uint >(ncells);

    // The index has one slot per state; it is rebuilt when the next run starts
    state_index_ready = false;

}

// Diseases and references ----------------------------------------------------

inline int Model::nh_global(const Track & t, int local) const
{
    int offset = 0;
    for (const auto & track : tracks)
    {
        if (track.id == t.id)
            break;
        offset += static_cast< int >(track.def->get_n_transitions());
    }
    return offset + local;
}

inline const Model::Track & Model::nh_local(int global, int & local) const
{
    local = global;
    if (local >= 0)
        for (const auto & track : tracks)
        {
            const int n = static_cast< int >(track.def->get_n_transitions());
            if (local < n)
                return track;
            local -= n;
        }

    throw std::range_error(
        "There is no transition " + std::to_string(global) + "."
    );
}

inline void Model::nh_disease_error(int disease, const char * fn) const
{

    if (tracks.empty())
        throw std::logic_error("The model has no natural history.");

    std::string names;
    for (const auto & n : get_diseases())
        names += (names.empty() ? "\"" : ", \"") + n + "\"";

    if (disease == all_diseases)
        throw std::logic_error(
            std::string(fn) + " needs the disease: the model has several (" +
            names + "). Pass its index (disease_id(name))."
        );

    throw std::range_error(
        std::string(fn) + ": there is no disease " + std::to_string(disease) +
        " (the model's diseases are " + names + ")."
    );

}

inline int Model::nh_resolve_name(
    const std::string & name,
    const std::string & what
) const
{

    if (!name.empty())
        return disease_id(name);

    if (tracks.size() == 1u)
        return 0;

    std::string names;
    for (const auto & n : get_diseases())
        names += (names.empty() ? "\"" : ", \"") + n + "\"";

    throw std::logic_error(
        what + " does not name its disease, and the model has several (" +
        names + "): add .of(\"<disease>\")."
    );

}

inline void Model::nh_resolve()
{

    if (tracks.empty())
        return;

    for (auto & r : nh_reactions)
    {
        if (r.on.kind != Trigger::Kind::entered)
            continue;
        r.track = nh_resolve_name(
            r.on.disease,
            "The reaction to entering stage " + std::to_string(r.on.target)
        );
        tracks[static_cast< size_t >(r.track)].def->check_stage(r.on.target, "trigger's");
    }

    // The tools must name their diseases when there are several; a layer
    // without targets carries every disease
    for (auto & t : tools)
    {
        t->targets_.resolve(*this, "The tool \"" + t->get_name() + "\"", true);
        db.reg.tool_mask[static_cast< size_t >(t->get_id())] = t->targets_.mask;
    }
    for (auto & l : layers)
        l->targets_.resolve(*this, "The layer \"" + l->get_name() + "\"", false);

    const size_t nt = st ? st->get_n_transitions() : 0u;
    nh_st_when.assign(nt, 0);
    nh_st_anchor.assign(nt, 0);
    for (size_t t = 0u; t < nt; ++t)
    {

        const auto & tr = st->get_transition(static_cast< int >(t));
        if (!tr.when.empty())
        {
            nh_st_when[t] = nh_resolve_name(
                tr.when_disease,
                "The status transition \"" + tr.name + "\" (its when())"
            );
            for (int s : tr.when)
                tracks[static_cast< size_t >(nh_st_when[t])].def->check_stage(s, "when");
        }

        if ((tr.kind == TransitionKind::after) && (tr.duration.get_anchor() >= 0))
        {
            nh_st_anchor[t] = nh_resolve_name(
                tr.duration.get_anchor_disease(),
                "The duration of the status transition \"" + tr.name +
                "\" (its since())"
            );
            tracks[static_cast< size_t >(nh_st_anchor[t])].def->check_stage(
                tr.duration.get_anchor(), "anchor"
            );
        }

    }

}

inline void DiseaseTargets::resolve(
    const Model & model,
    const std::string & owner,
    bool required
)
{

    if (names.empty())
    {
        if (required && (model.get_n_diseases() > 1u))
        {
            std::string all;
            for (const auto & n : model.get_diseases())
                all += (all.empty() ? "\"" : ", \"") + n + "\"";
            throw std::logic_error(
                owner + " does not name the diseases it acts on, and the model "
                "has several (" + all + "): add them with add_target()."
            );
        }
        mask = ~uint64_t(0);
        return;
    }

    mask = 0u;
    for (const auto & n : names)
    {
        if (!model.get_n_diseases())
            throw std::logic_error(owner + " targets \"" + n + "\", but the model has no disease.");
        mask |= uint64_t(1) << model.disease_id(n);
    }

}

// Accessors -------------------------------------------------------------------

inline std::vector< std::string > Model::get_diseases() const
{
    std::vector< std::string > names;
    for (const auto & track : tracks)
        names.push_back(track.nh->get_name());
    return names;
}

inline int Model::disease_id(std::string_view name) const
{

    for (size_t d = 0u; d < tracks.size(); ++d)
        if (tracks[d].nh->get_name() == name)
            return static_cast< int >(d);

    std::string names;
    for (const auto & n : get_diseases())
        names += (names.empty() ? "\"" : ", \"") + n + "\"";

    throw std::range_error(
        "The model has no disease \"" + std::string(name) + "\". Its diseases "
        "are: " + (names.empty() ? std::string("none") : names) + "."
    );

}

inline const std::vector< std::string > & Model::get_stages(int disease) const
{
    return nh_disease(disease, "get_stages()").nh->get_stages();
}

inline int Model::get_stage(const Agent & agent, int disease) const
{
    return nh_stage_of(nh_disease(disease, "get_stage()"), agent.state);
}

inline int Model::get_status(const Agent & agent) const
{
    if (tracks.empty())
        throw std::logic_error("The model has no natural history.");
    return nh_status_of(agent.state);
}

inline int Model::transition_id(std::string_view name, int disease) const
{

    if (tracks.empty())
        throw std::logic_error("The model has no natural history.");

    int found = -1;
    int n_found = 0;
    for (const auto & track : tracks)
    {
        if ((disease != all_diseases) && (track.id != disease))
            continue;
        for (size_t t = 0u; t < track.def->get_n_transitions(); ++t)
            if (track.def->get_transition(static_cast< int >(t)).name == name)
            {
                found = nh_global(track, static_cast< int >(t));
                ++n_found;
            }
    }

    if (n_found > 1)
        throw std::logic_error(
            "Several diseases have a transition named \"" + std::string(name) +
            "\": pass the disease (transition_id(name, disease_id(...)))."
        );

    if (n_found == 1)
        return found;

    if ((disease != all_diseases) && (static_cast< size_t >(disease) >= tracks.size()))
        nh_disease(disease, "transition_id()");

    if (st)
        return st->transition_id(name);

    throw std::logic_error(
        "The model has no transition \"" + std::string(name) + "\"."
    );

}

inline int Model::get_stage_entered(size_t id, int stage, int disease) const
{
    const Track & track = nh_disease(disease, "get_stage_entered()");
    track.nh->check_stage(stage, "");
    const int row = track.row.at(id);
    if (row < 0)
        return (stage == 0) ? 0 : -1;
    return track.rows[static_cast< size_t >(row) * (1u + track.nh->get_n_stages()) +
        1u + static_cast< size_t >(stage)];
}

inline int Model::get_stage_until(size_t id, int disease) const
{
    return nh_disease(disease, "get_stage_until()").until.at(id);
}

inline int Model::get_episode(size_t id, int disease) const
{
    const Track & track = nh_disease(disease, "get_episode()");
    const int row = track.row.at(id);
    if (row < 0)
        return 0;
    return track.rows[static_cast< size_t >(row) * (1u + track.nh->get_n_stages())];
}

inline int * Model::Track::record(size_t id)
{

    const size_t width = 1u + nh->get_n_stages();
    if (row[id] < 0)
    {
        // No infection yet, in the first stage since day 0
        row[id] = static_cast< int >(rows.size() / width);
        rows.resize(rows.size() + width, -1);
        rows[rows.size() - width]      = 0;
        rows[rows.size() - width + 1u] = 0;
    }

    return rows.data() + static_cast< size_t >(row[id]) * width;

}

inline bool Model::Track::removed(int stage) const
{
    // A stage that takes agents out of every layer and has no exit (the
    // built-in "Removed (other disease)" follows a removal elsewhere)
    if ((stage == removed_other) || (nh->get_contact(stage) > 0.0))
        return false;
    for (size_t t = 0u; t < nh->get_n_transitions(); ++t)
        if (nh->is_source(static_cast< int >(t), stage))
            return false;
    return true;
}

// Moving agents from code ------------------------------------------------------

inline bool Model::fire(int transition, Agent & agent)
{

    if (tracks.empty())
        throw std::logic_error("fire() needs a natural history.");

    if (nh_is_status_tr(transition))
    {

        const int k = nh_status_tr(transition);
        if (!st || (static_cast< size_t >(k) >= st->get_n_transitions()))
            throw std::range_error(
                "There is no status transition " + std::to_string(k) + "."
            );

        if (!st->is_source(k, nh_status_of(agent.state)))
            return false;

        for (const auto & track : tracks)
            if (track.removed(nh_stage_of(track, agent.state)))
                return false;

        nh_propose(tracks[0u], transition, agent);
        return true;

    }

    int local = 0;
    const Track & track = nh_local(transition, local);
    const auto & tr = track.nh->get_transition(local);

    if (tr.kind == TransitionKind::infection)
        throw std::logic_error(
            "An infection is not fired; seed the virus instead."
        );

    if (!track.nh->is_source(local, nh_stage_of(track, agent.state)))
        return false;

    // One infection at a time (several-diseases.md, Section 4.3)
    if (agent.virus && (agent.virus->disease != track.id))
    {
        bool carrier = (tr.to >= 0) && track.nh->is_carrier(tr.to);
        for (int to : tr.branch_to)
            carrier = carrier || track.nh->is_carrier(to);
        if (carrier)
            return false;
    }

    nh_propose(track, local, agent);
    return true;

}

inline void Model::react(Trigger on, int transition)
{
    if (tracks.empty())
        throw std::logic_error("react() needs a natural history.");
    nh_reactions.push_back({std::move(on), transition, nullptr, 0});
}

inline void Model::react(
    Trigger on,
    std::function< void(Agent &, Model &) > fun
)
{
    if (tracks.empty())
        throw std::logic_error("react() needs a natural history.");
    nh_reactions.push_back({std::move(on), -1, std::move(fun), 0});
}

inline void Model::set_initial(Assign a, int disease)
{
    const int d = nh_disease(disease, "set_initial()").id;
    tracks[static_cast< size_t >(d)].initial = std::move(a);
}

// The run ---------------------------------------------------------------------

inline void Model::nh_reset()
{

    const size_t n = population.size();

    bool first_timed = false;
    for (auto & track : tracks)
    {
        track.until.assign(n, INT_MAX);
        track.row.assign(n, -1);
        track.rows.clear();
        first_timed = first_timed || (track.nh->stages[0u].after >= 0);
    }
    nh_status_until.assign(st ? n : 0u, INT_MAX);
    nh_status_sod.assign(st ? n : 0u, 0);
    nh_pending = -1;
    nh_pending_track = 0;

    // Everyone starts in the first stage of each disease and the first
    // status, on day 0 (per agent: its stages' durations, then its status's)
    const int s_after = st ? st->stages[0u].after : -1;
    if (first_timed || (s_after >= 0))
        for (auto & p : population)
        {
            const size_t id = static_cast< size_t >(p.id);
            for (auto & track : tracks)
            {
                const int after = track.nh->stages[0u].after;
                if (after >= 0)
                    track.until[id] = track.nh->transitions[after].duration.draw(p, *this);
            }
            if (s_after >= 0)
                nh_status_until[id] = st->transitions[s_after].duration.draw(p, *this);
        }

    if (use_queuing && (nh_queue[0u] == 1))
        for (size_t id = 0u; id < n; ++id)
            queue.shift(id, 1);

}

inline void Model::nh_step(Agent & p)
{

    // An infection (one, across the diseases), then each disease's stage,
    // then the status (judged on the stages at the start of the day). With
    // one disease, nh_build_states() calls its stage update directly.
    const unsigned int c = p.state;
    const int infected = (p.virus == nullptr) ? nh_step_infection(p, c) : -1;
    for (auto & track : tracks)
        if (track.id != infected)
            nh_step_stage(track, p, c, false);

    if (st)
        nh_step_status(p, c);

}

inline int Model::nh_step_infection(Agent & p, unsigned int c)
{

    // Not out of the layers
    if (!state_excluded.empty() && state_excluded[c])
        return -1;

    // The diseases on the default sampler share one draw (one
    // trial, -1); each custom sampler is a trial of its own. They are tried
    // in a random order; the first infection ends the trial.
    uint64_t on_default = 0u;
    auto & order = nh_infection_order;
    order.clear();
    for (const auto & track : tracks)
    {
        const int inf = track.nh->stages[static_cast< size_t >(nh_stage_of(track, c))].infection;
        if (inf < 0)
            continue;
        if (track.nh->transitions[static_cast< size_t >(inf)].sampler)
            order.push_back(track.id);
        else
            on_default |= (uint64_t(1) << track.id);
    }

    if (on_default != 0u)
        order.push_back(-1);

    for (size_t i = order.size(); i > 1u; --i)
        std::swap(order[i - 1u], order[runif_index(static_cast< uint32_t >(i))]);

    for (const int trial : order)
    {

        Virus * v = nullptr;
        if (trial < 0)
        {
            nh_trial_mask = on_default;
            v = sample_virus(&p);
        }
        else
        {
            const Track & track = tracks[static_cast< size_t >(trial)];
            const int inf = track.nh->stages[static_cast< size_t >(nh_stage_of(track, c))].infection;
            nh_trial_mask = uint64_t(1) << trial;
            v = track.nh->transitions[static_cast< size_t >(inf)].sampler(p, *this);
            if ((v != nullptr) && (v->disease != trial))
            {
                nh_trial_mask = ~uint64_t(0);
                throw std::logic_error(
                    "The infection sampler of \"" + track.nh->get_name() +
                    "\" returned a virus of another disease (\"" +
                    v->get_name() + "\")."
                );
            }
        }
        nh_trial_mask = ~uint64_t(0);

        if (v != nullptr)
        {
            const Track & track = tracks[static_cast< size_t >(v->disease)];
            nh_propose(
                track,
                track.nh->stages[static_cast< size_t >(nh_stage_of(track, c))].infection,
                p, v
            );
            return v->disease;
        }

    }

    return -1;

}

inline bool Model::nh_step_stage(
    Track & track,
    Agent & p,
    unsigned int c,
    bool one
)
{

    const int s = nh_stage_of(track, c);
    const StageInfo & si = track.nh->stages[static_cast< size_t >(s)];

    // The values that apply to the agent: its virus's (a strain's), or the
    // track's (with one disease the virus is always of it; see nh_of())
    const auto & v = p.virus;
    const NaturalHistory & h =
        (v && v->nh_def && (one || (v->disease == track.id))) ?
        *v->nh_def : *track.nh;

    // 1. Infection by contact (one disease; not out of the layers)
    if (
        one && (si.infection >= 0) && (p.virus == nullptr) &&
        (state_excluded.empty() || !state_excluded[c])
    )
    {
        const auto & tr = h.transitions[si.infection];
        Virus * v = tr.sampler ? tr.sampler(p, *this) : sample_virus(&p);
        if (v != nullptr)
        {
            if (v->disease != track.id)
                throw std::logic_error(
                    "The infection sampler of \"" + track.nh->get_name() +
                    "\" returned a virus of another disease (\"" +
                    v->get_name() + "\")."
                );
            nh_propose(track, si.infection, p, v);
            return true;
        }
    }

    // 2. The rates compete: one draw picks one, or none
    const size_t nr = si.rates.size();
    if (nr > 0u)
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
            nh_propose(track, si.rates[static_cast< size_t >(which)], p);
            return true;
        }

    }

    // 3. The duration, if it ended (a deadline: the rates went first)
    if ((si.after >= 0) && (track.until[static_cast< size_t >(p.id)] <= today()))
    {
        nh_propose(track, si.after, p);
        return true;
    }

    return false;

}

inline void Model::nh_step_status(Agent & p, unsigned int c)
{

    const int k = nh_status_of(c);
    const StageInfo & ki = st->stages[static_cast< size_t >(k)];

    // A status transition restricted to some stages (of its disease; with
    // one disease, the stage is read once)
    const int s0 = nh_stage_of(tracks[0u], c);
    auto applies = [this, c, s0](int t) -> bool {
        const auto & tr = st->transitions[static_cast< size_t >(t)];
        if (tr.when.empty())
            return true;
        const int s = (tracks.size() == 1u) ? s0 :
            nh_stage_of(tracks[static_cast< size_t >(nh_st_when[static_cast< size_t >(t)])], c);
        return std::find(tr.when.begin(), tr.when.end(), s) != tr.when.end();
    };

    const size_t nkr = ki.rates.size();
    if (nkr > 0u)
    {

        reserve_tmp(nkr);
        for (size_t r = 0u; r < nkr; ++r)
        {
            const auto & tr = st->transitions[ki.rates[r]];
            array_double_tmp[r] = applies(ki.rates[r]) ? tr.value.eval(p, *this) : 0.0;
        }

        const int which = roulette(static_cast< epiworld_fast_uint >(nkr), this);
        if (which >= 0)
        {
            nh_propose(tracks[0u], -2 - ki.rates[static_cast< size_t >(which)], p);
            return;
        }

    }

    if (
        (ki.after >= 0) &&
        (nh_status_until[static_cast< size_t >(p.id)] <= today()) &&
        applies(ki.after)
    )
        nh_propose(tracks[0u], -2 - ki.after, p);

}

inline void Model::nh_propose(
    const Track & track,
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
            nh_of(track, p).transitions[transition]);

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

    // The event carries the target (a stage or a status) and its track;
    // events_run() makes it a cell when it applies it
    nh_pending = transition;
    nh_pending_track = track.id;

    const int s = nh_stage_of(track, p.state);
    if (status)
        p.change_state(*this, static_cast< epiworld_fast_uint >(to), Queue::NoOne);
    else if (tr.kind == TransitionKind::infection)
        p.set_virus(*this, *virus, to, Queue::NoOne);
    else if (track.nh->is_carrier(s) && !track.nh->is_carrier(to))
        p.rm_virus(*this, to, Queue::NoOne);
    else
        p.change_state(*this, static_cast< epiworld_fast_uint >(to), Queue::NoOne);

    nh_pending = -1;
    nh_pending_track = 0;

}

inline void Model::nh_enter(
    Track & track,
    Agent & p,
    unsigned int stage_old,
    int transition
)
{

    const size_t id = static_cast< size_t >(p.id);
    const size_t ns = track.nh->get_n_stages();
    const int s = nh_stage_of(track, p.state);

    // A new infection starts a new episode: its entry days start afresh
    const bool infected =
        (transition == NH_SEED) ||
        (!track.nh->is_carrier(static_cast< int >(stage_old)) && track.nh->is_carrier(s));

    int * record  = track.record(id);
    int * entered = record + 1;
    if (infected)
    {
        ++record[0];
        std::fill(entered, entered + ns, -1);
    }

    entered[s] = today();

    const int after = track.nh->stages[static_cast< size_t >(s)].after;
    if (after >= 0)
    {

        const auto & d = nh_of(track, p).transitions[after].duration;
        int base = today();
        if ((d.get_anchor() >= 0) && (entered[d.get_anchor()] >= 0))
            base = entered[d.get_anchor()];

        track.until[id] = base + d.draw(p, *this);

    }
    else
        track.until[id] = INT_MAX;

}

inline void Model::nh_react(Agent & p, int transition, int track, int stage_old)
{

    if (nh_reactions.empty())
        return;

    const size_t id = static_cast< size_t >(p.id);
    const Track & t = tracks[static_cast< size_t >(track)];
    const int s = nh_stage_of(t, p.state);

    // A stage transition's id in the model
    const int fired = (transition >= 0) ? nh_global(t, transition) : transition;

    // The status at the start of the day (the agent changed today)
    const int status_sod = st ? nh_status_sod[id] : 0;

    for (const auto & r : nh_reactions)
    {

        const bool hit = (r.on.kind == Trigger::Kind::entered) ?
            ((r.track == track) && (s == r.on.target) && (stage_old != s)) :
            (fired == r.on.target);

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

    for (auto & track : tracks)
    {

        if (!track.initial.is_set())
            continue;

        // The agents left without a virus, in the disease's first stage
        auto & buf = assign_buffers;
        buf.eligible.resize(population.size());
        for (auto & p : population)
            buf.eligible[static_cast< size_t >(p.id)] =
                (p.virus == nullptr) && (nh_stage_of(track, p.state) == 0);

        track.initial.draw(*this, buf);

        for (const auto & x : buf.out)
        {

            const size_t i = x.first;
            const int s    = x.second;

            if (
                (s < 0) || (static_cast< size_t >(s) >= track.def->get_n_stages()) ||
                track.nh->is_carrier(s)
            )
                throw std::logic_error(
                    "set_initial(): the label " + std::to_string(s) + " is not a "
                    "stage of \"" + track.nh->get_name() + "\" without the "
                    "virus (an initial stage is a stage, e.g., "
                    "assign::shares({{R, 0.1}}); carriers are seeded with the "
                    "virus)."
                );

            nh_pending = NH_INITIAL;
            nh_pending_track = track.id;
            population[i].change_state(
                *this, static_cast< epiworld_fast_uint >(s), Queue::NoOne
            );
            nh_pending = -1;
            nh_pending_track = 0;

        }

        events_run();

    }

}

inline epiworld_double Model::nh_infection_prob(
    Agent & target,
    const Virus & v,
    const Track & track
)
{
    const int tr = track.nh->stages[
        static_cast< size_t >(nh_stage_of(track, target.state))
    ].infection;
    if (tr < 0)
        return 0.0;

    // The infector's virus says how infectious it is
    const auto & h = v.nh_def ? *v.nh_def : *track.nh;
    return h.transitions[tr].value.eval(target, &v, *this);
}

inline std::string Model::draw_natural_history() const
{

    if (tracks.empty())
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

    // One disease: "stages" and "s<i>"; several: one subgraph each
    for (const auto & track : tracks)
    {
        const std::string d = (tracks.size() == 1u) ? "" : std::to_string(track.id);
        subgraph(
            "stages" + d, "Stages: " + track.nh->get_name(),
            track.nh->draw_mermaid((tracks.size() == 1u) ? "s" : "d" + d + "s")
        );
    }
    if (st)
        subgraph("statuses", "Statuses", st->draw_mermaid("q"));

    out +=
        "    classDef infectious fill:#f8d7da,stroke:#b02a37\n"
        "    classDef carrier fill:#fff3cd,stroke:#997404\n"
        "    classDef away fill:#e2e3e5,stroke:#6c757d\n";

    return out;

}

#endif
