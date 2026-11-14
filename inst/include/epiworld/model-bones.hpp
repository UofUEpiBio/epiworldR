#ifndef EPIWORLD_MODEL_BONES_HPP
#define EPIWORLD_MODEL_BONES_HPP

#include <memory>
#include <vector>
#include <functional>
#include <string>
#include <map>
#include "config.hpp"

#include "naturalhistory-bones.hpp"
#include "agent-bones.hpp"
#include "virus-bones.hpp"
#include "tool-bones.hpp"
#include "database-bones.hpp"
#include "queue-bones.hpp"
#include "globalevent-bones.hpp"
#include "contacttracing-bones.hpp"
#include "postsampling-bones.hpp"
#include "assign.hpp"

class AgentsSample;

class AdjList;

inline std::function<void(size_t,Model*)> make_save_run(
    std::string fmt = "%03lu-episimulation.csv",
    bool total_hist = true,
    bool virus_info = false,
    bool virus_hist = false,
    bool tool_info = false,
    bool tool_hist = false,
    bool transmission = false,
    bool transition = false,
    bool reproductive = false,
    bool generation = false,
    bool active_cases = false,
    bool outbreak_size = false,
    bool hospitalizations = false
    );

/**
 * @brief Read-only view of a contiguous run of agent ids.
 *
 * @details Returned by `Model::get_agents_in_state()`. It points into the
 * model's index, so it is invalidated by the next step (or anything else that
 * moves agents between states).
 */
class AgentIdsView {
private:
    const size_t * first = nullptr;
    size_t n = 0u;
public:
    AgentIdsView() = default;
    AgentIdsView(const size_t * first_, size_t n_) : first(first_), n(n_) {}
    const size_t * begin() const { return first; }
    const size_t * end() const { return first + n; }
    size_t size() const { return n; }
    bool empty() const { return n == 0u; }
    size_t operator[](size_t i) const { return first[i]; }
};

/**
 * @brief Core class of epiworld.
 *
 * The model class provides the wrapper that puts together `Agent`, `Virus`, and
 * `Tools`.
 *
 * virus and human sequence is represented as numeric vectors (if needed.)
 */
class Model {
    friend class Agent;
    friend class AgentsSample;
    friend class DataBase;
    friend class Queue;
    friend class layer::Network;
    friend class layer::Mixing;

protected:

    std::string name = ""; ///< Name of the model

    DataBase db = DataBase(*this);

    std::vector< Agent > population = {};

    /**
     * @brief Creates the population: `n` agents, in the baseline state.
     * @details The models' constructors call it; the population is set once,
     * when the model is built (a different size is a different model).
     */
    void set_population_size(size_t n);

    /**
     * @name Auxiliary variables for AgentsSample iterators
     *
     * @details These variables+objects are used by the AgentsSample
     * class for building efficient iterators over agents. The idea is to
     * reduce the memory allocation, so only during the first call of
     * AgentsSample::AgentsSample(Model) these vectors are allocated.
     */
    ///@{
    std::vector< Agent * > sampled_population;
    size_t sampled_population_n = 0u;
    std::vector< size_t > population_left;
    size_t population_left_n = 0u;
    ///@}

    /**
     * @name Agents features
     *
     * @details Optionally, a model can include an external data source
     * pointing to agents information. The data can then be access through
     * the `Agent::operator()` method.
     *
     */
    ///@{
    double * agents_data = nullptr;
    size_t agents_data_ncols = 0u;
    ///@}

    std::vector< VirusPtr > viruses = {};
    std::vector< ToolPtr > tools = {};
    std::vector< Assign > virus_assign = {}; ///< [virus] Who gets it at each reset
    std::vector< Assign > tool_assign = {};  ///< [tool] Who gets it at each reset
    AssignBuffers assign_buffers;            ///< The draws' work space; never copied

    std::shared_ptr< epi_xoshiro256ss > engine = std::make_shared< epi_xoshiro256ss >();

    epiworld_double runifd_a = 0.0;
    epiworld_double runifd_b = 1.0;
    std::normal_distribution<>       rnormd      =
        std::normal_distribution<>(0.0);
    std::gamma_distribution<>        rgammad     =
        std::gamma_distribution<>();
    std::lognormal_distribution<>    rlognormald =
        std::lognormal_distribution<>();
    std::exponential_distribution<>  rexpd       =
        std::exponential_distribution<>();
    std::binomial_distribution<> rbinomd         =
        std::binomial_distribution<>();
    int rbinomd_n = 1;
    epiworld_double rbinomd_fast_lambda = 0.0;
    bool rbinomd_use_poisson = false;
    std::negative_binomial_distribution<> rnbinomd =
        std::negative_binomial_distribution<>();
    std::geometric_distribution<> rgeomd          =
        std::geometric_distribution<>();
    std::poisson_distribution<> rpoissd           =
        std::poisson_distribution<>();

    // Parameter values, by position, and each name's position. Parameters are
    // never removed, so positions (ParamId) stay valid. The comparator is
    // transparent (std::less<>) so names can be looked up as
    // std::string_view without building a std::string.
    std::vector< epiworld_double > param_values;
    std::map<std::string, size_t, std::less<> > param_index;
    uint32_t param_layout_id = new_param_layout_id();
    epiworld_fast_uint ndays = 0;
    Progress pb;

    std::vector< UpdateFun >    state_fun = {};                  ///< Functions to update states
    std::vector< std::string >        states_labels = {};              ///< Labels of the states

    epiworld_fast_uint nstates = 0u;

    bool verbose     = true;
    int current_date = 0;

    // True while run() is driving the day loop, so get_ndays() is the run's
    // actual horizon. Not copied: a copy of the model is not in a run.
    bool running = false;

    void set_running(bool value); ///< Also tells the layers (see `ContactLayer`)

    void dist_tools();
    void dist_virus();

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed =
        std::chrono::duration<epiworld_double,std::micro>::zero();
    epiworld_fast_uint n_replicates = 0u;
    int last_seed = 0;
    void chrono_start();
    void chrono_end();

    std::vector<GlobalEventPtr> globalevents;

    /// Contact layers (see `add_layer()`), deep-copied with the model
    std::vector< std::shared_ptr< ContactLayer > > layers;

    /**
     * @name Network layers
     *
     * @details The network layers among `layers`, in order. Susceptible agents
     * pull from all of them, the queue and pushing follow all of their ties,
     * and the agents-by-state index counts an agent's ties across them (its
     * degree).
     */
    ///@{
    std::vector< layer::Network * > network_layers;
    std::vector< layer::Mixing * > group_layers; ///< Within-group mixing layers (the queue follows them)
    std::vector< size_t > layer_contacts; ///< Scratch of `sample_virus()`
    bool pull_only = false; ///< See `has_pull_only_layer()`
    bool queue_blind = false; ///< See `has_queue_blind_layer()`
    int sampled_layer = -1; ///< See `take_sampled_layer()`

    /// Adds the exposures of `contacts` to the arrays; returns their number.
    size_t add_exposures(
        Agent * agent,
        const ContactLayer * layer,
        const size_t * contacts,
        size_t n,
        size_t n_exposures
    );

    /// One draw over the first `n_exposures` exposures.
    Virus * draw_exposure(size_t n_exposures);
    void layers_attach(); ///< Points the layers at this model, lists the networks
    void layers_classify(); ///< Sets `pull_only`, `queue_blind` and `group_layers`
    void layers_detach(); ///< Leaves the layers without a model (they may outlive it)
    void add_layer_ptr(std::shared_ptr< ContactLayer > layer);
    size_t network_degree(size_t id) const; ///< Ties of agent `id`, all networks

    /**
     * Whether some layer's contacts are not seen from both ends or not by the
     * queue: a directed network (a tie is kept by its source only) or a layer
     * that is not a network. The model then pulls, and visits every agent
     * instead of using the queue.
     */
    bool has_pull_only_layer() const { return pull_only; }

    /**
     * Whether some layer's contacts are not seen by the queue (see
     * `ContactLayer::supports_queue()`): a directed network, or a mixing layer
     * with a contact matrix. The model then visits every agent.
     */
    bool has_queue_blind_layer() const { return queue_blind; }

    /// Agent `i` moved from `from` to `to` in a within-group mixing layer.
    void group_member_moved(
        size_t i,
        const std::vector< size_t > * from,
        const std::vector< size_t > * to
    );

    /// A network layer changed the ties of `i` and `j` by `di` and `dj`.
    void network_tie_changed(size_t i, int di, size_t j, int dj);

    /// A network layer swapped two ties (see `layer::Network::swap_ties()`).
    void network_ties_swapped(
        size_t i, size_t j, size_t k, size_t l, bool directed
    );
    ///@}

    Queue queue;
    bool use_queuing   = true;
    size_t sim_id = 0u;
    void set_sim_id(size_t id);

    PostSamplingFun post_sampling_fun = nullptr; ///< See set_post_sampling()
    bool post_sampling_on = false;                     ///< post_sampling_fun is set
    PostSamplingScratch post_sampling_scratch;         ///< Never copied
    void post_sampling_dispatch();
    void post_sampling_prepare_scratch(); ///< Sizes the scratch for the population

    /**
     * @name Agents by state
     *
     * @details All agent ids live in one array, `state_order`, grouped by
     * state: the agents in state `s` are `state_order[state_start[s]]` up to
     * (not including) `state_order[state_start[s + 1]]`, in no particular
     * order, and `state_member_pos[i]` is where agent `i` sits. One contiguous
     * block of N ids, whatever the number of states. Moving an agent from
     * state `a` to state `b` walks it across the blocks in between, one swap
     * per block boundary: O(|a - b|).
     *
     * Per state, the index also keeps the sum of the members' degrees, how
     * many of them carry a virus, and the sum of those carriers' degrees.
     *
     * Built in `reset()` and kept current in `events_run()` -- the only place
     * where an agent's state or virus changes during a run -- and by the
     * network layers' edits. The degree sums let the transmission step
     * compare the cost of pushing and pulling in O(number of states).
     */
    ///@{
    std::vector< size_t > state_order;       ///< [N] Agent ids, grouped by state
    std::vector< size_t > state_start;       ///< [nstates + 1] Where each state's block starts
    std::vector< size_t > state_member_pos;  ///< [agent] Position in state_order
    std::vector< unsigned int > agent_state;    ///< [agent] Copy of Agent::state, compact for scans
    std::vector< char > agent_carrier;          ///< [agent] 1 if Agent::virus is set, compact for scans
    std::vector< size_t > state_degree;
    std::vector< size_t > state_carriers;
    std::vector< size_t > state_carrier_degree;
    bool state_index_ready = false;

    void state_index_build();
    void state_index_update(Agent * p, unsigned int state_old, bool had_virus);
    void state_index_move(size_t id, unsigned int state_old, unsigned int state_new);
    AgentIdsView state_index_members(size_t state) const;
    void state_index_degree(size_t id, int delta);
    ///@}

    /**
     * @name Infectious and excluded states
     * @details Derived from the natural history and the statuses (see
     * `nh_build_states()`). The per-state masks are built at `reset()`.
     */
    ///@{
    std::vector< int > infectious_states;     ///< Empty: every state
    std::vector< int > excluded_states;
    std::vector< bool > state_not_infectious; ///< [state] Empty if every state is infectious
    std::vector< char > state_excluded;       ///< [state]
    void state_masks_build();
    ///@}

    /**
     * @name Network transmission
     *
     * @details See `set_transmission_mode()` and
     * `model-meat-transmission.hpp`. The push scratch space is per model and
     * is not copied: a copy sizes its own the first time it pushes.
     */
    ///@{
    TransmissionMode transmission_mode = TransmissionMode::automatic;
    TransmissionMode transmission_mode_last = TransmissionMode::pull;
    double transmission_kappa = EPI_DEFAULT_TRANSMISSION_KAPPA;

    struct PushTarget {
        size_t id;
        double odds;                // Sum of p / (1 - p) over the contacts
        unsigned int n_certain;     // Contacts with p >= 1
        Virus * candidate;    // Infector's virus drawn so far
        int layer;                  // Layer of the candidate's contact
    };

    std::vector< char > push_pushable;       ///< [state] Infected by a default infection only
    std::vector< char > push_source_ok;      ///< [state] Infectious, and some state is pushable
    std::vector< int >  push_slot;           ///< [agent] Index in push_targets, or -1
    std::vector< PushTarget > push_targets;
    std::vector< uint64_t > push_visit;      ///< [agent bit] To update after a push
    std::vector< uint64_t > push_sources;    ///< [agent bit] Carriers that push this step

    bool transmission_prepare();
    bool transmission_choose_push() const;
    void transmission_push();
    template<bool Record> void transmission_push_impl();
    void transmission_update_others();
    ///@}

    /**
     * @name Natural history
     *
     * @details See `naturalhistory-bones.hpp`. The first virus with a natural
     * history installs it (`add_virus()`); the model's copy of that virus
     * shares it. It is shared by the model's copies (it does not change once
     * installed); what it stores per agent is the model's, and is reset with
     * it.
     *
     * What the engine keeps for one disease is a **track** (see
     * `several-diseases.md`). A model has one track for now; the functions
     * below take theirs, so that a second disease is a second track.
     */
    ///@{
    struct Track {
        std::shared_ptr< const NaturalHistory > nh; ///< The disease
        std::vector< int > until;  ///< [agent] Day the stage's duration ends
        std::vector< int > row;    ///< [agent] Its record in `rows` (-1: none yet)
        /**
         * Per agent that has left the first stage (`row`), `1 + stages`
         * ints: its infections so far, then the day it entered each stage in
         * this infection (-1: not). An agent without a record has been in
         * the first stage since day 0, with no infection; the records are
         * made as agents move, so a run touches memory for the agents that
         * do.
         */
        std::vector< int > rows;
        std::vector< char > queue; ///< [cell] Queue mode: 0 none, 1 self, 2 everyone
        Assign initial;            ///< See set_initial()

        int * record(size_t id);   ///< The agent's record, made on first use
        bool removed(int stage) const; ///< No contacts and no exit
    };
    std::vector< Track > tracks;

    /// The track of a virus (its disease), and of an agent's infection
    /// (without one, the first)
    Track & nh_track(const Virus *) { return tracks[0u]; }
    const Track & nh_track(const Virus *) const { return tracks[0u]; }

    /// The agent's cell in a track. The first track's is kept in the agent.
    unsigned int nh_cell_of(const Track &, const Agent & p) const { return p.state; }

    int nh_pending = -1;            ///< Transition of the next event (see _add_event())

    /// Event: a seeded infection, and an initial placement (status transitions
    /// are `-2 - index`, see `TransitionRef`)
    static constexpr int NH_SEED = std::numeric_limits< int >::min();
    static constexpr int NH_INITIAL = std::numeric_limits< int >::min() + 1;

    /// The statuses (`nullptr`: one status, "Free")
    std::shared_ptr< const Statuses > st = nullptr;
    std::vector< int > nh_status_until; ///< [agent] Day the status's duration ends (with statuses)
    std::vector< int > nh_status_sod;   ///< [agent] Status at the start of the day it last changed (with statuses)

    struct Reaction {
        Trigger on;
        int transition;                          ///< Fired, or (-1) ...
        std::function< void(Agent &, Model &) > fun; ///< ... called
    };
    std::vector< Reaction > nh_reactions;

    size_t nh_nstatus() const { return st ? st->get_n_stages() : 1u; }
    // Without statuses a cell is its stage (no division on the hot path)
    int nh_stage_of(unsigned int cell) const
    {
        return static_cast< int >(st ? cell / st->get_n_stages() : cell);
    }
    int nh_status_of(unsigned int cell) const
    {
        return static_cast< int >(st ? cell % st->get_n_stages() : 0u);
    }
    unsigned int nh_cell(int stage, int status) const
    {
        return static_cast< unsigned int >(stage) * static_cast< unsigned int >(nh_nstatus()) +
            static_cast< unsigned int >(status);
    }
    static bool nh_is_status_tr(int t) { return (t <= -2) && (t > NH_INITIAL); }
    static int nh_status_tr(int t) { return -2 - t; }
    int nh_db_id(int t) const
    {
        return nh_is_status_tr(t) ?
            static_cast< int >(tracks[0u].nh->get_n_transitions()) + nh_status_tr(t) : t;
    }
    void nh_build_states();
    void nh_react(Agent & p, int transition, int stage_old);

    void nh_install(std::shared_ptr< NaturalHistory > def);
    void nh_reset();
    void nh_step(Agent & p);
    bool nh_step_stage(Track & track, Agent & p, unsigned int cell); ///< nh_step()'s stage part
    void nh_step_status(Agent & p, unsigned int cell, int stage); ///< nh_step()'s status part

    /// The natural history whose values apply to `p` in a track: its
    /// virus's (a strain's values), or (no virus) the track's
    const NaturalHistory & nh_of(const Track & track, const Agent & p) const
    {
        const auto & v = p.get_virus();
        return (v && v->nh_def) ? *v->nh_def : *track.nh;
    }
    void nh_propose(int transition, Agent & p, Virus * virus = nullptr);
    void nh_enter(Track & track, Agent & p, unsigned int state_old, int transition);
    void nh_queue_move(Agent & p, unsigned int state_old);
    void nh_initial_place();
    epiworld_double nh_infection_prob(Agent & target, const Virus & v);
    ///@}

    /**
     * @brief Variables used to keep track of the events
     * to be made regarding viruses.
     */
    std::vector< Event > events = {};
    epiworld_fast_uint nactions = 0u;

    /**
     * @brief Construct a new Event object
     *
     * @param agent_ Agent over which the action will be called
     * @param virus_ Virus pointer included in the action
     * @param tool_ Tool pointer included in the action
     * @param new_state_ New state of the agent
     * @param queue_ Change in the queue
      * @param action_ Action to execute when processing the event
     */
    void _add_event(
        Agent * agent_,
        VirusPtr virus_,
        ToolPtr tool_,
        epiworld_fast_int new_state_,
        epiworld_fast_int queue_,
        EventAction action_,
        int layer_ = -1
        );

    /**
     * @name Default event handlers
     *
     * These member functions implement the default behavior for each
     * event action (add/remove virus or tool, or change state).
     */
    ///@{
    void _event_add_virus(Event & a);
    void _event_add_tool(Event & a);
    void _event_rm_virus(Event & a);
    void _event_rm_tool(Event & a);
    void _event_change_state(Event & a);
    ///@}

    /**
     * @name Tool Mixers
     *
     * These functions combine the effects tools have to deliver
     * a single effect. For example, wearing a mask, been vaccinated,
     * and the immune system combine together to jointly reduce
     * the susceptibility for a given virus.
     *
     */
    virtual epiworld_double susceptibility_reduction_mixer(
        Agent * agent, VirusPtr & virus
    );
    virtual epiworld_double transmission_reduction_mixer(
        Agent * agent, VirusPtr & virus
    );
    virtual epiworld_double recovery_enhancer_mixer(
        Agent * agent, VirusPtr & virus
    );
    virtual epiworld_double death_reduction_mixer(
        Agent * agent, VirusPtr & virus
    );

    /**
     * @brief Advanced usage: Makes a copy of data and returns it as undeleted pointer
     *
     * @param copy
     */
    virtual std::unique_ptr<Model> clone_ptr();

    /// For the models built on this class, which set their population.
    Model();

public:

    /**
     * @name Scratch arrays for sampling exposures
     * @details Update functions store each exposure's probability in
     * `array_double_tmp` and its virus in `array_virus_tmp`, then draw with
     * `roulette()`, which uses two doubles per exposure. An agent can have any
     * number of exposures across its layers, so call `reserve_tmp()` before
     * storing them.
     */
    ///@{
    std::vector< epiworld_double > array_double_tmp =
        std::vector< epiworld_double >(1024u * 2u);
    std::vector< Virus * > array_virus_tmp =
        std::vector< Virus * >(1024u);
    std::vector< int > array_layer_tmp =
        std::vector< int >(1024u); ///< Layer of each exposure

    /// Makes room for `n` exposures in the scratch arrays.
    void reserve_tmp(size_t n);
    ///@}

    /// A model with `n` agents.
    explicit Model(size_t n);
    Model(const Model & m);
    Model(Model && m);
    Model & operator=(const Model & m);

    virtual ~Model() { layers_detach(); };

    DataBase & get_db();
    const DataBase & get_db() const;
    epiworld_double operator()(std::string_view pname) const;

    size_t size() const;

    /**
     * @name Random number generation
     *
     * @param eng Random number generator
     * @param s Seed
     */
    ///@{
    void set_rand_engine(std::shared_ptr< epi_xoshiro256ss > & eng);
    std::shared_ptr< epi_xoshiro256ss > & get_rand_endgine();
    void seed(size_t s);
    void set_rand_norm(epiworld_double mean, epiworld_double sd);
    void set_rand_unif(epiworld_double a, epiworld_double b);
    void set_rand_exp(epiworld_double lambda);
    void set_rand_gamma(epiworld_double alpha, epiworld_double beta);
    void set_rand_lognormal(epiworld_double mean, epiworld_double shape);
    void set_rand_binom(int n, epiworld_double p);
    void set_rand_nbinom(int n, epiworld_double p);
    void set_rand_geom(epiworld_double p);
    void set_rand_poiss(epiworld_double lambda);
    epiworld_double rnorm();
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma();
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double rexp();
    epiworld_double rexp(epiworld_double lambda);
    epiworld_double rlognormal();
    epiworld_double rlognormal(epiworld_double mean, epiworld_double shape);

    /**
     * @brief Draw from the currently configured uniform distribution.
     * @return A random draw from the configured uniform distribution.
     * @details
     * These uniform draws make use of Lemire's algorithm for fast
     * uniform integer generation, which is both faster and more accurate than
     * the common `std::uniform_int_distribution` approach.
     * 
     * 
     * Lemire, D. (2019). Fast Random Integer Generation in an Interval.
     * ACM Trans. Model. Comput. Simul., 29(1), 3:1-3:12.
     * <https://doi.org/10.1145/3230636>
     */
    epiworld_double runif();
    epiworld_double runif(epiworld_double a, epiworld_double b);
    int runif_int(int a, int b);
    uint32_t runif_index(uint32_t n);
    /**
     * @brief Draw from the currently configured binomial distribution.
     * @details When `EPI_FAST_BINOM` is enabled (default), this uses
     * `rpoiss(lambda)` with `lambda = n * p` after `set_rand_binom(n, p)` in the
     * rare-event regime `p <= 0.01` and `n * p * p <= 0.1`. Define
     * `EPI_NO_FAST_BINOM` before including epiworld to disable this behavior.
     * @return A random draw from the configured binomial distribution, or from
     * its Poisson approximation when the fast path is active.
     */
    int rbinom();
    /**
     * @brief Draw from a binomial distribution with parameters `n` and `p`.
     * @details When `EPI_FAST_BINOM` is enabled (default), this uses
     * `rpoiss(lambda)` with `lambda = n * p` in the rare-event regime
     * `p <= 0.01` and `n * p * p <= 0.1`. This preserves the mean and is often
     * substantially faster in practice while remaining very accurate in that
     * region. Define `EPI_NO_FAST_BINOM` before including epiworld to disable
     * this behavior and force the exact binomial draw.
     * @param n Number of trials.
     * @param p Success probability.
     * @return A random draw from the binomial distribution, or from its
     * Poisson approximation when the fast path is active.
     */
    int rbinom(int n, epiworld_double p);
    int rnbinom();
    int rnbinom(int n, epiworld_double p);
    int rgeom();
    int rgeom(epiworld_double p);
    int rpoiss();
    int rpoiss(epiworld_double lambda);
    ///@}

    /**
     * @brief Sample from a set of probabilities stored in array_double_tmp.
     * @details Uses a cumulative probability approach: draws a uniform random
     * number and walks through array_double_tmp[0..n-1], accumulating
     * probabilities until the draw is exceeded. If no event fires, returns n
     * (meaning "none of the above").
     * @param n Number of probability entries in array_double_tmp to consider.
     * @return Index in [0, n] of the sampled event (n = no event).
     */
    size_t sample_from_probs(size_t n);

    /**
     * @name Add Virus/Tool to the model
     *
     * This is done before the model has been initialized.
     *
     * @param v Virus to be added
     * @param t Tool to be added
     * @param preval Initial prevalence (initial state.) It can be
     * specified as a proportion (between zero and one,) or an integer
     * indicating number of individuals.
     */
    ///@{
    /**
     * @brief Adds a virus, and who gets it at the start of each run.
     * @details `a` gives each agent without a virus whose stage can be
     * infected a label: the stage it starts in (a carrier stage), or, for the
     * forms without labels (e.g., `assign::share()`), the infection's target.
     * By default, no one gets it (see `assign.hpp`).
     *
     * Each call records a new variant, the founder of a lineage, whose id is
     * the number of viruses added before it (adding the same virus twice
     * records two; a removed virus keeps its record, so ids do not move).
     * Mutations during a run add variants after them, which the next run
     * (or the next `add_virus()`) drops.
     */
    void add_virus(const Virus & v, Assign a = Assign());

    /// Adds a tool, and who gets it at the start of each run (default: no one).
    void add_tool(const Tool & t, Assign a = Assign());

    /// Replaces who gets virus `id` (e.g., the virus of a built-in model).
    void assign_virus(size_t id, Assign a);

    /// Replaces who gets tool `id`.
    void assign_tool(size_t id, Assign a);

    void rm_virus(size_t virus_pos);
    void rm_tool(size_t tool_pos);
    ///@}

    /**
     * @name Accessing population of the model
     */
    ///@{
    std::vector< Agent > & get_agents(); ///< Returns a reference to the vector of agents.

    Agent & get_agent(size_t i);

    /**
     * @brief Each agent's stage in a disease, and each agent's status
     * (indices into `get_stages()` and the statuses).
     * @param disease The disease (`disease_id()`); `all_diseases`, the
     * default, is the one disease of a model with one.
     */
    ///@{
    std::vector< int > get_agents_stages(int disease = all_diseases) const;
    std::vector< int > get_agents_statuses() const;
    ///@}




    /**
     * @name Functions to run the model
     *
     * @param seed Seed to be used for Pseudo-RNG.
     * @param ndays Number of days (steps) of the simulation.
     * @param fun In the case of `run_multiple`, a function that is called
     * after each experiment.
     *
     */
    ///@{
    void update_state();
    void mutate_virus();
    virtual void next();
    virtual Model & run(
        epiworld_fast_uint ndays,
        int seed = -1
    ); ///< Runs the simulation (after initialization)
    Model & run_multiple( ///< Multiple runs of the simulation
        epiworld_fast_uint ndays,
        epiworld_fast_uint nexperiments,
        int seed_ = -1,
        std::function<void(size_t,Model*)> fun = make_save_run(),
        bool reset = true,
        bool verbose = true,
        int nthreads = 1
        );
    ///@}

    size_t get_n_viruses() const; ///< Number of viruses in the model
    size_t get_n_tools() const; ///< Number of tools in the model
    epiworld_fast_uint get_ndays() const;

    /**
     * @brief True while `run()` is driving the day loop.
     *
     * @details It is set from just before `reset()` until the last day is
     * done, so global events (including their `reset()`) can rely on
     * `get_ndays()` being the run's horizon -- even when it is zero. It is
     * false when the day loop is driven by hand (`reset()`, then the steps
     * called directly), where `get_ndays()` means nothing.
     */
    bool is_running() const;
    epiworld_fast_uint get_n_replicates() const;
    size_t get_sim_id() const;
    void set_ndays(epiworld_fast_uint ndays);
    bool get_verbose() const;
    Model & verbose_off();
    Model & verbose_on();
    int today() const; ///< The current time of the model


    /**
     * @brief Wrapper of `DataBase::write_data`
     *
     * @param fn_virus_info Filename. Information about the virus.
     * @param fn_virus_hist Filename. History of the virus.
     * @param fn_tool_info Filename. Information about the tool.
     * @param fn_tool_hist Filename. History of the tool.
     * @param fn_total_hist   Filename. Aggregated history (state)
     * @param fn_transmission Filename. Transmission history.
     * @param fn_transition   Filename. Markov transition history.
     * @param fn_reproductive_number Filename. Case by case reproductive number
     * @param fn_generation_time Filename. Generation time data.
     * @param fn_active_cases Filename. Active cases data.
     * @param fn_outbreak_size Filename. Outbreak size data.
     * @param fn_hospitalizations Filename. Hospitalization data.
     */
    void write_data(
        std::string fn_virus_info,
        std::string fn_virus_hist,
        std::string fn_tool_info,
        std::string fn_tool_hist,
        std::string fn_total_hist,
        std::string fn_transmission,
        std::string fn_transition,
        std::string fn_reproductive_number,
        std::string fn_generation_time,
        std::string fn_active_cases,
        std::string fn_outbreak_size,
        std::string fn_hospitalizations
        ) const;


    /**
     * @brief Copy of the model parameters, by name.
     * @details To change a parameter, use `set_param()`.
     */
    std::map<std::string, epiworld_double> params() const;

    /**
     * @brief Reset the model
     *
     * @details Resetting the model will:
     * - clear the database
     * - reset the contact layers (see `ContactLayer`)
     * - re-distribute tools
     * - re-distribute viruses
     * - set the date to 0
     *
     */
    virtual void reset();
    const Model & print(bool lite = false) const;

    /**
     * @brief The number of the engine's states: one per pair of stage (of the
     * virus's natural history) and status (of the model's statuses), stage
     * major (`state = stage * n_statuses + status`). Use `get_stage()` and
     * `get_status()` for an agent's.
     */
    size_t get_n_states() const;

    /**
     * @name Natural history
     * @brief The disease's stages and the declared transitions between them.
     */
    ///@{

    /// The natural history of the model's viruses (`nullptr` if none has one).
    const NaturalHistory * get_natural_history() const
    {
        return tracks.empty() ? nullptr : tracks[0u].nh.get();
    }

    /**
     * @name The model's diseases
     * @details A disease is a natural history; its strains share it (see
     * `several-diseases.md`). A model has one for now.
     */
    ///@{
    size_t get_n_diseases() const { return tracks.size(); }
    std::vector< std::string > get_diseases() const; ///< Their names, by index
    int disease_id(std::string_view name) const;     ///< Throws if there is none
    ///@}

    /**
     * @brief Fires a declared transition on an agent, applied at the end of
     * the step (or of the global event).
     * @details Nothing happens if, when it is applied, the agent is not in one
     * of the transition's sources. An infection is not fired: seed instead.
     * @return Whether the agent was in a source when fired.
     */
    bool fire(int transition, Agent & agent);

    /**
     * @brief Places agents left without a virus in stages at the start of
     * each run (e.g., `assign::shares({{R, 0.1}})`: 10% Recovered).
     * @details The labels are stages without the virus; the eligible agents
     * are those left without a virus, in the first stage. Replaces what was
     * set before.
     */
    void set_initial(Assign a);

    /// Day agent `id` entered `stage` in its current (or last) infection (-1: not).
    int get_stage_entered(size_t id, int stage) const;

    /// Day agent `id`'s stage duration ends (`INT_MAX`: none).
    int get_stage_until(size_t id) const;

    /// Infections agent `id` has had in the current run.
    int get_episode(size_t id) const;

    /// The stage and the status of an agent.
    int get_stage(const Agent & agent) const;
    int get_status(const Agent & agent) const;

    /**
     * @brief Sets the model's statuses (see `Statuses`): where agents are,
     * and the transitions between them. Between runs only.
     */
    void set_statuses(Statuses s);

    /// The statuses (`nullptr` if the model has one, "Free").
    const Statuses * get_statuses() const { return st.get(); }

    /**
     * @brief When `on` happens to an agent, fire `transition` on it (or call
     * `fun`), in the same step.
     * @details Reactions link the axes: e.g., entering Rash fires
     * "Quarantined -> Isolated". Within a day, stage changes apply first, then
     * status changes, then reactions (which see the resulting status; `.was()`
     * reads the one at the start of the day).
     */
    void react(Trigger on, int transition);
    void react(Trigger on, std::function< void(Agent &, Model &) > fun);

    /// A transition's id (stage or status), by name.
    int transition_id(std::string_view name) const;
    ///@}

    /**
     * @brief Ids of the agents currently in a state.
     *
     * @details The list is kept up to date as the model runs, so looking up
     * who is in a state costs nothing (no scan of the population). The ids are
     * in no particular order. The index is built when a run starts, so this is
     * available once `run()` (or `run_multiple()`) has been called, and it
     * reflects the model at its current step.
     *
     * @param state The state code.
     * @return A view of the ids (iterable, with `size()` and `[]`); it is
     * invalidated by the next step.
     * @throws std::logic_error if the model has not been run yet (or its
     * population changed since).
     * @throws std::range_error if `state` is not a state of the model.
     */
    AgentIdsView get_agents_in_state(epiworld_fast_uint state) const;

    /**
     * @brief One flag per agent: 1 if the agent carries a virus.
     * @details A compact copy of "has a virus", one byte per agent, so a scan
     * of an agent's neighbors can rule most of them out without loading their
     * `Agent` (a cache miss each). It is kept current from the moment a run
     * starts. Returns `nullptr` before that (and after the population
     * changes): callers then have to look at the agents.
     */
    const char * get_agent_carrier() const
    {
        return state_index_ready ? agent_carrier.data() : nullptr;
    }

    /**
     * @brief Draws the virus, if any, that infects `agent` from its contacts.
     * @details This is where contacts become an infection. Every contact that
     * carries a virus is an exposure with probability
     * `(1 - susceptibility reduction of agent) * probability of infecting *
     * (1 - transmission reduction of the contact)`, and one draw
     * (`roulette()`) picks at most one infector. A contact listed twice is
     * two exposures.
     * With a `layer`, the probabilities follow its weight and the viruses it
     * admits (see `ContactLayer`), and the layer is recorded as the one of the
     * transmission.
     * @param agent Agent without a virus.
     * @param contacts Ids of the contacts.
     * @param n Number of contacts.
     * @param layer The layer the contacts came from (none: weight 1, any
     * virus, no layer recorded).
     * @return The virus transmitted, or `nullptr` if none is.
     */
    Virus * sample_virus(
        Agent * agent,
        const size_t * contacts,
        size_t n,
        const ContactLayer * layer = nullptr
    );

    /**
     * @brief The layer of the last transmission drawn by `sample_virus()`
     * (-1 if none), which the next `Agent::set_virus()` records and clears.
     */
    int take_sampled_layer()
    {
        const int layer = sampled_layer;
        sampled_layer = -1;
        return layer;
    }

    /**
     * @brief Draws the virus, if any, that `agent` catches from its contacts
     * in the model's layers.
     * @details The contacts that carry a virus, layer by layer in the order
     * the layers were added and each layer in its own order, are all
     * exposures of one draw (`sample_virus(agent, contacts, n, layer)` adds
     * them). Only contacts in an infectious stage (`stage::infectious`) are
     * sources. This is what a default infection does; a pushing step
     * (networks only) gives the same distribution.
     * @throws std::logic_error if `agent` carries a virus.
     */
    Virus * sample_virus(Agent * agent);

    /**
     * @name Network transmission mode
     *
     * @details Stages whose only drawn exit is a default infection (no
     * custom sampler) can be updated by pulling (each
     * susceptible agent scans its neighbors) or by pushing (each agent with a
     * virus adds its infection odds to its susceptible neighbors). Both give
     * the same distribution of who gets infected, and by whom; only the random
     * number stream differs. See `TransmissionMode`.
     *
     * With `"auto"` (the default) the model estimates the cost of each step,
     * as the ties to scan plus a fixed price per agent visited
     * (`EPI_TRANSMISSION_AGENT_COST` ties), for the carriers that can transmit
     * (push) and for the susceptible agents (pull). It pushes whenever the
     * push cost is no more than `kappa` times the pull cost, and pulls
     * otherwise; `kappa` only matters in this mode. The choice depends only on
     * the model's state, never on the queueing system, so turning queuing on or
     * off leaves results unchanged. A `kappa` below 1 makes pushing harder to
     * choose (it favors pulling); the default, 0.5, is below 1 because the
     * queue spares a pull the susceptibles with no infectious neighbor, which
     * the costs do not see. It was fitted together with the per-agent price
     * with `examples/20-transmission-benchmark`.
     *
     * Directed networks, and states with other update functions, always pull.
     * Set `"pull"` to reproduce the random streams of epiworld <= 0.15.
     *
     * @param mode `"auto"`, `"push"`, or `"pull"` (or the enum).
     * @param kappa Relative cost threshold used by `"auto"`: a finite,
     * non-negative number (default `EPI_DEFAULT_TRANSMISSION_KAPPA`, 0.5).
     * @throws std::invalid_argument for an unknown mode.
     * @throws std::range_error for a negative or infinite `kappa`.
     */
    ///@{
    Model & set_transmission_mode(
        TransmissionMode mode,
        double kappa = EPI_DEFAULT_TRANSMISSION_KAPPA
    );
    Model & set_transmission_mode(
        std::string_view mode,
        double kappa = EPI_DEFAULT_TRANSMISSION_KAPPA
    );
    TransmissionMode get_transmission_mode() const;
    /// The mode used in the most recent step (`push` or `pull`).
    TransmissionMode get_last_transmission_mode() const;
    /// The threshold used by `"auto"`.
    double get_transmission_kappa() const;

    /**
     * @brief The sums `"auto"` compares to choose between push and pull.
     * @details Totals over the states that can be sources of infection (the
     * carriers) and over the states updated by the default susceptible sampler
     * (the susceptibles), as they stand now. All zero before the model has run.
     */
    struct TransmissionSums {
        size_t carrier_degree;      ///< Sum of the carriers' degrees
        size_t susceptible_degree;  ///< Sum of the susceptibles' degrees
        size_t carriers;            ///< Number of carriers that can transmit
        size_t susceptibles;        ///< Number of susceptible agents
    };
    TransmissionSums get_transmission_sums() const;
    ///@}


    /**
     * @name Setting and accessing parameters from the model
     *
     * @details Parameters are stored in a vector, in the order they were
     * added, with a map from each name to its position. Lookups by name take
     * the name as a `std::string_view` (no allocation) and search the map once.
     * Unknown names throw.
     *
     * Code that reads parameters for many agents every day can skip the
     * search:
     *
     * - `get_param_id()` returns the position (`ParamId`) of a parameter, and
     *   `par_at()` / `set_param_at()` read and write by position. Positions
     *   are valid for the model and all its copies.
     * - `ParamRef` (and the `EPI_PAR(model, "name")` macro) resolves a name
     *   once per model layout and caches the position, so it is safe to use
     *   with any model.
     *
     * The `par()` function members are aliases for `get_param()`.
     *
     * In the case of the function `read_params`, users can pass a file
     * listing parameters to be included in the model. Each line in the
     * file should have the following structure:
     *
     * ```
     * [name of parameter 1]: [value in double]
     * [name of parameter 2]: [value in double]
     * ...
     * ```
     *
     * The only condition for parameter names is that these do not include
     * a colon.
     *
     *
     * @param initial_val
     * @param pname Name of the parameter to add or to fetch
     * @param fn Path to the file containing parameters
     * @return The current value of the parameter
     * in the model.
     *
     */
    ///@{
    epiworld_double add_param(
        epiworld_double initial_val, std::string pname, bool overwrite = false
    );
    Model & read_params(std::string fn, bool overwrite = false);
    epiworld_double get_param(std::string_view pname) const;
    bool has_param(std::string_view pname) const;
    void set_param(std::string_view pname, epiworld_double val);
    epiworld_double par(std::string_view pname) const;
    ParamId get_param_id(std::string_view pname) const;
    /// Read on every step: the bounds check stays cheap (it throws from
    /// par_out_of_range()), so this can be inlined.
    epiworld_double par_at(ParamId id) const
    {
        if (id.idx >= param_values.size())
            par_out_of_range(id);
        return param_values[id.idx];
    }
    [[noreturn]] void par_out_of_range(ParamId id) const;
    void set_param_at(ParamId id, epiworld_double val);
    size_t get_n_params() const { return param_values.size(); };
    uint32_t get_param_layout_id() const { return param_layout_id; };
    ///@}

    void get_elapsed(
        std::string unit = "auto",
        epiworld_double * last_elapsed = nullptr,
        epiworld_double * total_elapsed = nullptr,
        std::string * unit_abbr = nullptr,
        bool print = true
    ) const;

    /**
     * @name Set the user data object
     *
     * @param names string vector with the names of the variables.
     */
    ///[@
    void set_user_data(std::vector< std::string > names);
    void add_user_data(epiworld_fast_uint j, epiworld_double x);
    void add_user_data(std::vector< epiworld_double > x);
    UserData & get_user_data();
    ///@}

    /**
     * @brief Set a global action
     *
     * @param fun A function to be called on the prescribed date
     * @param name Name of the action.
     * @param date Integer indicating when the function is called (see details)
     *
     * @details When date is less than zero, then the function is called
     * at the end of every day. Otherwise, the function will be called only
     * at the end of the indicated date.
     */
    void add_globalevent(
        std::function<void(Model*)> fun,
        std::string name = "A global action",
        int date = -99
        );

    void add_globalevent(
        GlobalEvent & action
    );

    GlobalEvent & get_globalevent(std::string name); ///< Retrieve a global action by name
    GlobalEvent & get_globalevent(size_t i); ///< Retrieve a global action by index
    bool has_globalevent(std::string_view name) const; ///< Whether a global action by that name exists
    size_t get_n_globalevents() const; ///< Number of global actions registered

    void rm_globalevent(std::string name); ///< Remove a global action by name
    void rm_globalevent(size_t i); ///< Remove a global action by index

    void run_globalevents();

    /**
     * @name Contact layers
     * @brief Sources of contacts of the model (see `ContactLayer`).
     * @{
     */

    /**
     * @brief Adds a contact layer to the model.
     * @details The model shares the layer with the caller (an edit through
     * `layer` is an edit of the model's layer), resets it at the start of
     * `reset()`, updates it at the end of every step, and clones it when the
     * model is copied (see `ContactLayer`). Layers are added between runs only
     * (a layer that should act only part of the time is added up front and
     * edited, e.g., by a global event).
     * @param layer The layer (e.g., from `layer::network()`). Its name must not
     * be in use by another layer, and it must not belong to another model.
     * @return A reference to the layer, as its own type.
     */
    template<typename TLayer>
    TLayer & add_layer(std::shared_ptr< TLayer > layer)
    {
        TLayer * ptr = layer.get();
        add_layer_ptr(std::move(layer));
        return *ptr;
    }

    size_t get_n_layers() const; ///< Number of contact layers
    bool has_layer(std::string_view name) const; ///< Whether a layer by that name exists
    size_t get_layer_id(std::string_view name) const; ///< Id of a layer by name (throws if absent)
    ContactLayer & get_layer(size_t id); ///< Retrieve a layer by id
    ContactLayer & get_layer(std::string_view name); ///< Retrieve a layer by name

    /// Retrieve a layer as its own type (unchecked cast; ids are cheaper than names)
    template<typename TLayer>
    TLayer & get_layer_as(size_t id)
    {
        return static_cast< TLayer & >(get_layer(id));
    }

    template<typename TLayer>
    TLayer & get_layer_as(std::string_view name)
    {
        return static_cast< TLayer & >(get_layer(name));
    }

    /// Calls `update()` on every layer (done by `next()` at the end of each step)
    void update_layers();

    /// The network layers of the model, in the order they were added.
    const std::vector< layer::Network * > & get_network_layers() const
    {
        return network_layers;
    }

    /** @} */

    void clear_state_set();

    /**
     * @name Queuing system
     * @details When queueing is on, the model will keep track of which agents
     * are either in risk of exposure or exposed. This then is used at each
     * step to act only on the aforementioned agents.
     *
     */
    ////@{
    void queuing_on(); ///< Activates the queuing system (default.)
    Model & queuing_off(); ///< Deactivates the queuing system.
    bool is_queuing_on() const; ///< Query if the queuing system is on.
    Queue & get_queue(); ///< Retrieve the `Queue` object.
    ///@}

    /**
     * @name Post-sampling callback
     * @details A callback that receives, for each infectious agent, the agents
     * it was in contact with during the step (see `PostSamplingFun`). The
     * models whose samplers report contacts (network pull and push, the
     * mixing models) collect them only while a callback is installed;
     * otherwise the cost is one check per sampling operation (per step in a
     * push). Pull and push report the same contacts. The callback is kept by
     * copies of the model, including the ones `run_multiple()` makes.
     *
     * Installing a callback replaces the previous one. Contact tracing does
     * not use it: a layer keeps its own log (`ContactLayer::enable_tracing()`).
     */
    ///@{
    Model & set_post_sampling(PostSamplingFun fun); ///< Install a callback (an empty function clears it).
    Model & clear_post_sampling(); ///< Remove the callback.
    bool has_post_sampling() const { return post_sampling_on; } ///< Is a callback installed?

    /// Reports that `infectious_id` was in contact with `contacted_id`. Samplers
    /// call it, only if `has_post_sampling()`.
    void register_sampled_contact(size_t infectious_id, size_t contacted_id);

    /// Same, for several infectious agents in contact with `contacted_id`.
    void register_sampled_contacts(
        const size_t * infectious_ids, size_t n, size_t contacted_id
    );
    ///@}

    const std::vector< VirusPtr > & get_viruses() const;
    const std::vector< ToolPtr > & get_tools() const;
    Virus & get_virus(size_t id);
    Virus & get_virus(std::string_view name);
    Tool & get_tool(size_t id);
    Tool & get_tool(std::string_view name);

    bool has_virus(std::string_view name) const;
    bool has_tool(std::string_view name) const;

    /**
     * @brief Set the agents data object
     *
     * @details The data should be an array with the data stored in a
     * column major order, i.e., by column.
     *
     * @param data_ Pointer to the first element of an array of size
     * `size() * ncols_`.
     * @param ncols_ Number of features included in the data.
     *
     */
    void set_agents_data(double * data_, size_t ncols_);
    double * get_agents_data();
    size_t get_agents_data_ncols() const;

    /**
     * @brief Set the name object
     *
     * @param name
     */
    void set_name(std::string name);
    std::string get_name() const;

    bool operator==(const Model & other) const;
    bool operator!=(const Model & other) const {return !operator==(other);};

    /**
     * @brief Executes the stored action
     *
     * @param model_ Model over which it will be executed.
     */
    void events_run();

    /**
     * @brief Draws a mermaid diagram of the model.
     * @param model The model to draw.
     * @param fn_output The name of the file to write the diagram.
     * If empty, the diagram will be printed to the standard output.
     * @param self Whether to allow self-transitions.
     */
    void draw(
        DiagramType diagram_type = DiagramType::Mermaid,
        const std::string & fn_output = "",
        bool self = false
    );

    /**
     * @brief The model's natural history as a mermaid flowchart: the stages
     * and, if the model has statuses, the statuses, each in a subgraph (see
     * `NaturalHistory::draw_mermaid()`). Drawn from the declarations, so it
     * needs no run.
     * @throws std::logic_error if the model has no natural history.
     */
    std::string draw_natural_history() const;

    /**
     * @brief Record a hospitalization event for an agent.
     * 
     * @param agent Reference to the agent being hospitalized.
     * 
     * @details
     * This is a wrapper for `DataBase::record_hospitalization()`.
     * For each hospitalization, the method records:
     * - The current date from the model
     * - The virus ID from the agent's virus
     * - For each tool the agent has, a separate record with weight = 1/N
     *   where N is the number of tools
     * - If the agent has no tools, a single record with tool_id = -1 and
     *   weight = 1.0
     */
    void record_hospitalization(Agent & agent);

    /**
     * @brief Get the full time series of hospitalization data.
     * 
     * @param date Output vector for dates.
     * @param virus_id Output vector for virus IDs.
     * @param tool_id Output vector for tool IDs.
     * @param count Output vector for counts (number of hospitalized individuals).
     * @param weight Output vector for summed weights (fractional contribution
     *   based on tool distribution).
     * 
     * @details
     * This is a wrapper for `DataBase::get_hospitalizations()`.
     * Returns the full time series of hospitalization data. For each unique 
     * (virus_id, tool_id) combination observed, returns an entry for every 
     * day from 0 to ndays-1.
     * 
     * The `count` vector contains the actual number of individuals hospitalized 
     * for that (date, virus_id, tool_id) combination. This is useful for 
     * answering questions like "how many total people were hospitalized?" 
     * regardless of their tools.
     * 
     * The `weight` vector contains fractional contributions: if an agent has N 
     * tools, each tool gets weight = 1/N. Summing weights across all tool_ids 
     * for a given date and virus_id gives the total number of hospitalizations.
     */
    void get_hospitalizations(
        std::vector<int> & date,
        std::vector<int> & virus_id,
        std::vector<int> & tool_id,
        std::vector<int> & count,
        std::vector<double> & weight
    ) const;


};

#endif
