#ifndef EPIWORLD_ASSIGN_HPP
#define EPIWORLD_ASSIGN_HPP

/**
 * @file assign.hpp
 * @brief How agents get something at the start of each run: a virus (and the
 * stage it puts them in), a tool, an initial stage, or a group of a mixing
 * layer. One vocabulary for all of them.
 *
 * @details An `Assign` gives each eligible agent a **label**, or none (-1),
 * drawn with the model's random numbers when a run starts:
 *
 * | Form | Labels |
 * |---|---|
 * | `assign::share(p)`, `assign::count(n)` | A random share (or number) of the eligible agents |
 * | `assign::shares({{l, p}, ...})`, `assign::counts({{l, n}, ...})` | Several labels: each takes its share (or number) of the eligible agents, in order, from those left |
 * | `assign::group_shares(layer, {p_g})`, `assign::group_counts(layer, {n_g})` | A share (or number) of each group of a mixing layer |
 * | `assign::agents(ids)` | The given agents |
 * | `assign::labels(v)` | `v[i]` for agent `i` (-1: none); the same every run |
 * | `assign::by(fun)` | `fun(agent, model)` for each eligible agent |
 *
 * What a label means, and who is eligible, depends on the use:
 *
 * | Use | Label | Eligible |
 * |---|---|---|
 * | `Model::add_virus(v, a)` | The stage the agent starts in (a carrier stage); a form without labels uses the infection's target | Agents without a virus whose stage can be infected |
 * | `Model::add_tool(t, a)` | -- (has the tool) | Everyone |
 * | `Model::set_initial(a)` | A stage without the virus | Agents left without a virus, in the first stage |
 * | `layer::Mixing::set_groups(a, n)` | A group, in `[0, n)` | Everyone |
 *
 * The placements are made in this order at every reset: the layers (so a
 * placement can use their groups), the viruses, the tools, the initial
 * stages, then the global events' set-up. The random forms are drawn again in
 * each run; `agents()` and `labels()` are fixed.
 */

class Assign;

namespace assign {
    /// The label of the forms that give none: the use's default (see above)
    constexpr int ANY = -2;
}

/// A function that gives an agent's label (-1: none).
using AssignFun = std::function< int(Agent &, Model &) >;

/**
 * @brief The work space of a draw. The model keeps one and reuses it at every
 * reset, so placing agents allocates nothing once its vectors have grown.
 * @details Copying it gives an empty one.
 */
struct AssignBuffers {

    /// In: one entry per agent, whether it can be given a label. A draw may
    /// clear entries (of the agents it labels).
    std::vector< char > eligible;

    /// Out: the agents given a label, with their label, in ascending id order.
    std::vector< std::pair< size_t, int > > out;

    /// The forms' work space (agent ids; the draws index it with 32 bits).
    std::vector< uint32_t > pool;

    AssignBuffers() = default;
    AssignBuffers(const AssignBuffers &) {}
    AssignBuffers & operator=(const AssignBuffers &) { return *this; }

};

/**
 * @brief Gives each eligible agent a label, or none (see `assign.hpp`).
 * @details Default-constructed, it gives no one anything.
 */
class Assign {
public:
    /// Appends the agents given a label, with their label, to `buf.out`,
    /// given who is eligible (`buf.eligible`).
    using DrawFun = std::function< void(Model &, AssignBuffers &) >;

private:
    DrawFun fun = nullptr;
    std::shared_ptr< const std::vector< int > > fixed_labels = nullptr;

public:
    Assign() = default;
    explicit Assign(
        DrawFun f,
        std::shared_ptr< const std::vector< int > > fixed_ = nullptr
    ) : fun(std::move(f)), fixed_labels(std::move(fixed_)) {}

    /// The labels, if they are the same every run (`assign::labels()`).
    const std::vector< int > * get_fixed() const { return fixed_labels.get(); }

    /// Whether it gives anyone anything.
    bool is_set() const { return static_cast< bool >(fun); }

    /**
     * @brief Draws the labels of the agents of `model`.
     * @param buf `buf.eligible` has one entry per agent: whether it can be
     * given a label (the draw may clear entries). On return, `buf.out` has the
     * agents given a label, with their label, in ascending id order.
     * @throws std::logic_error if the form names an agent that is not
     * eligible.
     */
    void draw(Model & model, AssignBuffers & buf) const;
};

#endif
