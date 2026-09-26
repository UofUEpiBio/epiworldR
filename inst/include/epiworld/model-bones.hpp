#ifndef EPIWORLD_MODEL_BONES_HPP
#define EPIWORLD_MODEL_BONES_HPP

#include <memory>
#include <vector>
#include <functional>
#include <string>
#include <map>
#include "config.hpp"

#include "agent-bones.hpp"
#include "virus-bones.hpp"
#include "viruses-bones.hpp"
#include "tool-bones.hpp"
#include "database-bones.hpp"
#include "queue-bones.hpp"
#include "globalevent-bones.hpp"
#include "contacttracing-bones.hpp"

template<typename TSeq>
class AgentsSample;

class AdjList;

template<typename TSeq = EPI_DEFAULT_TSEQ>
inline std::function<void(size_t,Model<TSeq>*)> make_save_run(
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

// template<typename TSeq>
// class VirusPtr;

// template<typename TSeq>
// class ToolPtr;

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
 * @tparam TSeq Type of sequence. In principle, users can build models in which
 * virus and human sequence is represented as numeric vectors (if needed.)
 */
template<typename TSeq>
class Model {
    friend class Agent<TSeq>;
    friend class AgentsSample<TSeq>;
    friend class DataBase<TSeq>;
    friend class Queue<TSeq>;

protected:

    std::string name = ""; ///< Name of the model

    DataBase<TSeq> db = DataBase<TSeq>(*this);

    std::vector< Agent<TSeq> > population = {};

    /// @brief Validates the arguments of `add_edge()` / `rm_edge()`.
    void check_edge_endpoints(size_t i, size_t j) const;

    bool using_backup = true;
    std::vector< Agent<TSeq> > population_backup = {};

    /**
     * @name Auxiliary variables for AgentsSample<TSeq> iterators
     *
     * @details These variables+objects are used by the AgentsSample<TSeq>
     * class for building efficient iterators over agents. The idea is to
     * reduce the memory allocation, so only during the first call of
     * AgentsSample<TSeq>::AgentsSample(Model<TSeq>) these vectors are allocated.
     */
    ///@{
    std::vector< Agent<TSeq> * > sampled_population;
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

    bool directed = false;

    std::vector< VirusPtr<TSeq> > viruses = {};
    std::vector< ToolPtr<TSeq> > tools = {};

    std::vector< Entity<TSeq> > entities = {};

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

    std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> rewire_fun;
    epiworld_double rewire_prop = 0.0;

    std::map<std::string, epiworld_double > parameters;
    epiworld_fast_uint ndays = 0;
    Progress pb;

    std::vector< UpdateFun<TSeq> >    state_fun = {};                  ///< Functions to update states
    std::vector< std::string >        states_labels = {};              ///< Labels of the states

    /** Function to distribute states. Goes along with the function  */
    std::function<void(Model<TSeq>*)> initial_states_fun = [](Model<TSeq> * /**/)
    -> void {};

    epiworld_fast_uint nstates = 0u;

    bool verbose     = true;
    int current_date = 0;

    // True while run() is driving the day loop, so get_ndays() is the run's
    // actual horizon. Not copied: a copy of the model is not in a run.
    bool running = false;

    void dist_tools();
    void dist_virus();
    void dist_entities();

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed =
        std::chrono::duration<epiworld_double,std::micro>::zero();
    epiworld_fast_uint n_replicates = 0u;
    int last_seed = 0;
    void chrono_start();
    void chrono_end();

    std::vector<GlobalEventPtr<TSeq>> globalevents;

    Queue<TSeq> queue;
    bool use_queuing   = true;
    size_t sim_id = 0u;
    void set_sim_id(size_t id);

    std::unique_ptr<ContactTracing> contact_tracing;
    bool use_contact_tracing = false;
    size_t contact_tracing_max_contacts = EPI_MAX_TRACKING;

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
     * where an agent's state or virus changes during a run -- and in
     * `add_edge()`/`rm_edge()`. The degree sums let the transmission step
     * compare the cost of pushing and pulling in O(number of states).
     */
    ///@{
    std::vector< size_t > state_order;       ///< [N] Agent ids, grouped by state
    std::vector< size_t > state_start;       ///< [nstates + 1] Where each state's block starts
    std::vector< size_t > state_member_pos;  ///< [agent] Position in state_order
    std::vector< unsigned int > agent_state;    ///< [agent] Copy of Agent::state, compact for scans
    std::vector< size_t > state_degree;
    std::vector< size_t > state_carriers;
    std::vector< size_t > state_carrier_degree;
    bool state_index_ready = false;

    void state_index_build();
    void state_index_update(Agent<TSeq> * p, unsigned int state_old, bool had_virus);
    void state_index_move(size_t id, unsigned int state_old, unsigned int state_new);
    AgentIdsView state_index_members(size_t state) const;
    void state_index_degree(Agent<TSeq> & p, size_t n_neighbors_before);
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
        Virus<TSeq> * candidate;    // Infector's virus drawn so far
    };

    std::vector< char > push_pushable;       ///< [state] Uses the default susceptible sampler
    std::vector< char > push_default;        ///< [state] ...and it is default_update_susceptible
    std::vector< char > push_excluded;       ///< [target state * nstates + source state]
    std::vector< char > push_source_ok;      ///< [state] Some pushable state accepts it as source
    std::vector< int >  push_slot;           ///< [agent] Index in push_targets, or -1
    std::vector< PushTarget > push_targets;
    std::vector< uint64_t > push_visit;      ///< [agent bit] To update after a push
    std::vector< uint64_t > push_sources;    ///< [agent bit] Carriers that push this step

    bool transmission_prepare();
    bool transmission_choose_push() const;
    void transmission_push();
    void transmission_update_others();
    ///@}

    /**
     * @brief Variables used to keep track of the events
     * to be made regarding viruses.
     */
    std::vector< Event<TSeq> > events = {};
    epiworld_fast_uint nactions = 0u;

    /**
     * @brief Construct a new Event object
     *
     * @param agent_ Agent over which the action will be called
     * @param virus_ Virus pointer included in the action
     * @param tool_ Tool pointer included in the action
     * @param entity_ Entity pointer included in the action
     * @param new_state_ New state of the agent
     * @param queue_ Change in the queue
      * @param action_ Action to execute when processing the event
     */
    void _add_event(
        Agent<TSeq> * agent_,
        VirusPtr<TSeq> virus_,
        ToolPtr<TSeq> tool_,
        Entity<TSeq> * entity_,
        epiworld_fast_int new_state_,
        epiworld_fast_int queue_,
          EventAction action_
        );

    /**
     * @name Default event handlers
     *
     * These member functions implement the default behavior for each
     * event action (add/remove virus, tool, entity, or change state).
     */
    ///@{
    void _event_add_virus(Event<TSeq> & a);
    void _event_add_tool(Event<TSeq> & a);
    void _event_add_entity(Event<TSeq> & a);
    void _event_rm_virus(Event<TSeq> & a);
    void _event_rm_tool(Event<TSeq> & a);
    void _event_rm_entity(Event<TSeq> & a);
    void _event_change_state(Event<TSeq> & a);
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
        Agent<TSeq> * agent, VirusPtr<TSeq> & virus
    );
    virtual epiworld_double transmission_reduction_mixer(
        Agent<TSeq> * agent, VirusPtr<TSeq> & virus
    );
    virtual epiworld_double recovery_enhancer_mixer(
        Agent<TSeq> * agent, VirusPtr<TSeq> & virus
    );
    virtual epiworld_double death_reduction_mixer(
        Agent<TSeq> * agent, VirusPtr<TSeq> & virus
    );

    /**
     * @brief Advanced usage: Makes a copy of data and returns it as undeleted pointer
     *
     * @param copy
     */
    virtual std::unique_ptr<Model<TSeq>> clone_ptr();

public:

    std::array<epiworld_double, 1024u * 2u> array_double_tmp;
    std::array<Virus<TSeq> *, 1024u * 2u> array_virus_tmp;

    Model();
    Model(const Model<TSeq> & m);
    Model(Model<TSeq> && m);
    Model<TSeq> & operator=(const Model<TSeq> & m);

    virtual ~Model() {};

    /**
     * @name Set the backup object
     * @details `backup` can be used to restore the entire object
     * after a run. This can be useful if the user wishes to have
     * individuals start with the same network from the beginning.
     * Building a new network (`agents_from_edgelist()` and friends, or
     * `agents_empty_graph()`) drops the backup.
     *
     */
    ///@{
    void set_backup();
    // void restore_backup();
    ///@}

    DataBase<TSeq> & get_db();
    const DataBase<TSeq> & get_db() const;
    epiworld_double operator()(std::string pname);

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
    void add_virus(Virus<TSeq> & v);
    void add_tool(Tool<TSeq> & t);
    void add_entity(Entity<TSeq> e);
    void rm_virus(size_t virus_pos);
    void rm_tool(size_t tool_pos);
    void rm_entity(size_t entity_id);
    ///@}

    /**
     * @brief Associate agents-entities from a file
     *
     * The structure of the file should be two columns separated by
     * space. The first column indexing between 0 and nagents-1, and the
     * second column between 0 and nentities - 1.
     *
     * @param fn Path to the file.
     * @param skip How many rows to skip.
     */
    void load_agents_entities_ties(std::string fn, int skip);

    /**
     * @brief Associate agents-entities from data
    */
    void load_agents_entities_ties(
        const std::vector<int> & agents_ids,
        const std::vector<int> & entities_ids
        );

    void load_agents_entities_ties(
        const int * agents_id,
        const int * entities_id,
        size_t n
        );

    /**
     * @name Accessing population of the model
     *
     * @details In a directed network (`directed = true`), a tie
     * `source -> target` is kept by its source only: `target` is one of
     * `source`'s neighbors, but not the other way around. An undirected tie is
     * kept at both ends. Update functions look at an agent's own neighbors (a
     * susceptible agent catches a virus from them), so the tie exposes the
     * source to the target. The target can infect the source, but not the
     * reverse. To have `i` infect `j`, give the tie as `j -> i`.
     *
     * A directed network always pulls (see `set_transmission_mode()`), and each
     * step updates (and offers mutations to) every agent, as with
     * `queuing_off()`. The queue flags the
     * neighbors of an agent that becomes infectious, and along a directed tie
     * those are not the agents it can infect. `write_edgelist()` returns the
     * ties as given, and `add_edge()`/`rm_edge()` refuse to edit a directed
     * network.
     *
     * @param fn std::string Filename of the edgelist file.
     * @param skip int Number of lines to skip in `fn`.
     * @param directed bool Whether the graph is directed or not.
     * @param size Size of the network.
     * @param al AdjList to read into the model. The network is directed if
     * `al` is.
     */
    ///@{
    void agents_from_adjlist(
        std::string fn,
        int size,
        int skip = 0,
        bool directed = false
        );

    void agents_from_edgelist(
        const std::vector< int > & source,
        const std::vector< int > & target,
        int size,
        bool directed
    );

    void agents_from_adjlist(AdjList al);

    bool is_directed() const; ///< Whether the network was built directed.

    std::vector< Agent<TSeq> > & get_agents(); ///< Returns a reference to the vector of agents.

    Agent<TSeq> & get_agent(size_t i);

    std::vector< epiworld_fast_uint > get_agents_states() const; ///< Returns a vector with the states of the agents.

    std::vector< Viruses_const<TSeq> > get_agents_viruses() const; ///< Returns a const vector with the viruses of the agents.

    std::vector< Viruses<TSeq> > get_agents_viruses(); ///< Returns a vector with the viruses of the agents.

    std::vector< Entity<TSeq> > & get_entities();

    Entity<TSeq> & get_entity(size_t entity_id, int * entity_pos = nullptr);
    const Entity<TSeq> & get_entity(size_t entity_id, int * entity_pos = nullptr) const;

    Model<TSeq> & agents_smallworld(
        epiworld_fast_uint n = 1000,
        epiworld_fast_uint k = 5,
        bool d = false,
        epiworld_double p = .01
        );
    /// Replaces the network with `n` agents and no ties (undirected).
    void agents_empty_graph(epiworld_fast_uint n = 1000);

    /**
     * @name Change the contact network of a model that already has one
     *
     * @details Unlike `agents_from_edgelist()` and friends, which build the
     * network up front, these edit it in place and may be called at any point,
     * including in the middle of a run from a global event -- a policy that
     * temporarily merges households, for instance. They keep the queueing system
     * in step as they go (see `Queue::notify_edge_added`), which is what makes
     * mid-run edits safe: the queue counts each agent's active neighbors, and a
     * tie appearing or disappearing underneath it would otherwise corrupt that
     * count and silently drop agents out of `update_state()`.
     *
     * Ties are undirected and are always changed at both ends. Existing
     * neighbors keep their relative order, so an edit never changes which
     * transmitter is sampled among the ties it left alone.
     *
     * @param i,j Ids of the two agents.
     * @throws std::range_error if an id is out of range.
     * @throws std::logic_error if `i == j`, or if the model is directed (these
     *         operate on both ends of a tie, which is meaningless there).
     */
    ///@{
    /**
     * @return `true` if the tie was created, `false` if the two were already
     *         tied. An intervention that has to withdraw its own ties later
     *         should record only the ones this returned `true` for, so it never
     *         removes a tie the model already had.
     */
    bool add_edge(size_t i, size_t j);

    bool rm_edge(size_t i, size_t j);  ///< @return `true` if a tie was removed.
    /// Whether `i` and `j` are tied. In a directed network, whether `i -> j` is.
    bool has_edge(size_t i, size_t j) const;
    ///@}

    /**
     * @brief Initialize agents using a Stochastic Block Model (SBM).
     *
     * Creates agents and connects them according to an SBM defined by
     * `block_sizes` and `mixing_matrix`.
     *
     * @param block_sizes Number of agents per block.
     * @param mixing_matrix Mixing matrix of size K*K; row sums give average
     *   expected degree per group.
     * @param row_major If `true`, matrix is row-major; otherwise column-major.
     * @return Reference to this Model.
     *
     * @see rgraph_sbm
     */
    Model<TSeq> & agents_sbm(
        const std::vector< size_t > & block_sizes,
        const std::vector< double > & mixing_matrix,
        bool row_major = true
        );
    ///@}

    /**
     * @name Initialize agents using a Bernoulli random graph
     * @param n Number of agents.
     * @param p Probability of tie formation.
     * @param d Whether the graph is directed or not.
     * @return Reference to this Model.
     */
    Model<TSeq> & agents_bernoulli(
        epiworld_fast_uint n,
        epiworld_double p,
        bool d = false
    );

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
    virtual Model<TSeq> & run(
        epiworld_fast_uint ndays,
        int seed = -1
    ); ///< Runs the simulation (after initialization)
    Model<TSeq> & run_multiple( ///< Multiple runs of the simulation
        epiworld_fast_uint ndays,
        epiworld_fast_uint nexperiments,
        int seed_ = -1,
        std::function<void(size_t,Model<TSeq>*)> fun = make_save_run<TSeq>(),
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
    size_t get_n_entities() const;
    void set_ndays(epiworld_fast_uint ndays);
    bool get_verbose() const;
    Model<TSeq> & verbose_off();
    Model<TSeq> & verbose_on();
    int today() const; ///< The current time of the model

    /**
     * @name Rewire the network preserving the degree sequence.
     *
     * @details In an undirected network, if {(i,j), (k,l)} -> {(i,l), (k,j)},
     * the reciprocal is also true, i.e., {(j,i), (l,k)} -> {(j,k), (l,i)}. In a
     * directed network only the sources' ties move, which keeps every agent's
     * in- and out-degree.
     *
     * The rewiring function runs on every step of a run, so it should change
     * ties only through `Agent::swap_neighbors()` (what `rewire_degseq()`
     * uses), `add_edge()`, or `rm_edge()`: these keep the queueing system in
     * step with the network.
     *
     * @param proportion Proportion of ties to be rewired.
     *
     * @result A rewired version of the network.
     */
    ///@{
    void set_rewire_fun(std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> fun);
    void set_rewire_prop(epiworld_double prop);
    epiworld_double get_rewire_prop() const;
    /// @brief Whether a rewiring function is set. `rewire()` calls it on every
    /// step, whatever the proportion.
    bool has_rewire_fun() const;
    void rewire();
    ///@}

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
     * @name Export the network data in edgelist form
     *
     * @param fn std::string. File name.
     * @param source Integer vector
     * @param target Integer vector
     *
     * @details When passing the source and target, the function will
     * write the edgelist on those.
     */
    ///@{
    void write_edgelist(
        std::string fn
        ) const;

    void write_edgelist(
        std::vector< int > & source,
        std::vector< int > & target
        ) const;
    ///@}

    std::map<std::string, epiworld_double> & params();

    /**
     * @brief Reset the model
     *
     * @details Resetting the model will:
     * - clear the database
     * - restore the population (if `set_backup()` was called before)
     * - re-distribute tools
     * - re-distribute viruses
     * - set the date to 0
     *
     */
    virtual void reset();
    const Model<TSeq> & print(bool lite = false) const;

    /**
     * @name Manage state (states) in the model
     *
     * @details
     *
     * The functions `get_state` return the current values for the
     * states included in the model. The function `set_state_function`
     * replaces the update function associated with an existing state.
     *
     * @param lab `std::string` Name of the state.
     *
     * @return `add_state*` returns the ID (index) of the registered state.
     * @return `set_state_function` returns a reference to the model.
     * @return `get_state_*` returns a vector of pairs with the
     * states and their labels.
     */
    ///@{
    epiworld_fast_int state_of(std::string_view name);
    epiworld_fast_int add_state(std::string lab, UpdateFun<TSeq> fun = nullptr);
    Model<TSeq> & set_state_function(epiworld_fast_uint state, UpdateFun<TSeq> fun = nullptr);
    Model<TSeq> & set_state_function(std::string_view name, UpdateFun<TSeq> fun = nullptr);
    const std::vector< std::string > & get_states() const;
    size_t get_n_states() const;
    const std::vector< UpdateFun<TSeq> > & get_state_fun() const;
    void print_state_codes() const;
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
     * @name Network transmission mode
     *
     * @details States whose update function is `default_update_susceptible`
     * or `sampler::make_update_susceptible()` can be updated by pulling (each
     * susceptible agent scans its neighbors) or by pushing (each agent with a
     * virus adds its infection odds to its susceptible neighbors). Both give
     * the same distribution of who gets infected, and by whom; only the random
     * number stream differs. See `TransmissionMode`.
     *
     * With `"auto"` (the default) the model pushes whenever the carriers'
     * ties are no more than `kappa` times the susceptibles' ties, and pulls
     * otherwise; `kappa` only matters in this mode. The choice depends only on the model's state, never on the
     * queueing system, so turning queuing on or off leaves results unchanged.
     * Because the queue already spares a pull the susceptibles with no
     * infectious neighbor -- which the rule does not see -- the default
     * `kappa` is 0.25 rather than 1 (tuned with
     * `examples/20-transmission-benchmark`).
     *
     * Directed networks, and states with other update functions, always pull.
     * Set `"pull"` to reproduce the random streams of epiworld <= 0.15.
     *
     * @param mode `"auto"`, `"push"`, or `"pull"` (or the enum).
     * @param kappa Relative cost threshold used by `"auto"`: a finite,
     * non-negative number (default `EPI_DEFAULT_TRANSMISSION_KAPPA`, 0.25).
     * @throws std::invalid_argument for an unknown mode.
     * @throws std::range_error for a negative or infinite `kappa`.
     */
    ///@{
    Model<TSeq> & set_transmission_mode(
        TransmissionMode mode,
        double kappa = EPI_DEFAULT_TRANSMISSION_KAPPA
    );
    Model<TSeq> & set_transmission_mode(
        std::string_view mode,
        double kappa = EPI_DEFAULT_TRANSMISSION_KAPPA
    );
    TransmissionMode get_transmission_mode() const;
    /// The mode used in the most recent step (`push` or `pull`).
    TransmissionMode get_last_transmission_mode() const;
    /// The threshold used by `"auto"`.
    double get_transmission_kappa() const;
    ///@}

    /**
     * @name Initial states
     *
     * @details These functions are called before the simulation starts.
     *
     * @param proportions_ Vector of proportions for each state.
     * @param queue_ Vector of queue for each state.
     */
    virtual Model<TSeq> & initial_states(
        std::vector< double > /*proportions_*/,
        std::vector< int > /*queue_*/
    ) {return *this;};

    /**
     * @name Setting and accessing parameters from the model
     *
     * @details Tools can incorporate parameters included in the model.
     * Internally, parameters in the tool are stored as pointers to
     * an std::map<> of parameters in the model. Using the `epiworld_fast_uint`
     * method directly fetches the parameters in the order these were
     * added to the tool. Accessing parameters via the `std::string` method
     * involves searching the parameter directly in the std::map<> member
     * of the model (so it is not recommended.)
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
    Model<TSeq> & read_params(std::string fn, bool overwrite = false);
    epiworld_double get_param(std::string pname);
    bool has_param(std::string_view pname) const;
    void set_param(std::string pname, epiworld_double val);
    epiworld_double par(std::string pname) const;
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
    UserData<TSeq> & get_user_data();
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
        std::function<void(Model<TSeq>*)> fun,
        std::string name = "A global action",
        int date = -99
        );

    void add_globalevent(
        GlobalEvent<TSeq> & action
    );

    GlobalEvent<TSeq> & get_globalevent(std::string name); ///< Retrieve a global action by name
    GlobalEvent<TSeq> & get_globalevent(size_t i); ///< Retrieve a global action by index
    bool has_globalevent(std::string_view name) const; ///< Whether a global action by that name exists
    size_t get_n_globalevents() const; ///< Number of global actions registered

    void rm_globalevent(std::string name); ///< Remove a global action by name
    void rm_globalevent(size_t i); ///< Remove a global action by index

    void run_globalevents();

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
    Model<TSeq> & queuing_off(); ///< Deactivates the queuing system.
    bool is_queuing_on() const; ///< Query if the queuing system is on.
    Queue<TSeq> & get_queue(); ///< Retrieve the `Queue` object.
    ///@}

    /**
     * @name Contact tracing
     * @details When contact tracing is on, the model will track contacts
     * between agents. Users must actively record contacts in their update
     * functions by calling `get_contact_tracing().add_contact(...)`.
     * Contact tracing is off by default.
     *
     * @param max_contacts Maximum number of contacts to track per agent
     * (default: EPI_MAX_TRACKING). Only used when turning tracing on.
     */
    ///@{
    Model<TSeq> & contact_tracing_on(size_t max_contacts = EPI_MAX_TRACKING); ///< Activates contact tracing.
    Model<TSeq> & contact_tracing_off(); ///< Deactivates contact tracing.
    bool is_contact_tracing_on() const; ///< Query if contact tracing is on.
    ContactTracing & get_contact_tracing(); ///< Retrieve the `ContactTracing` object.
    ///@}

    const std::vector< VirusPtr<TSeq> > & get_viruses() const;
    const std::vector< ToolPtr<TSeq> > & get_tools() const;
    Virus<TSeq> & get_virus(size_t id);
    Virus<TSeq> & get_virus(std::string_view name);
    Tool<TSeq> & get_tool(size_t id);
    Tool<TSeq> & get_tool(std::string_view name);

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

    bool operator==(const Model<TSeq> & other) const;
    bool operator!=(const Model<TSeq> & other) const {return !operator==(other);};

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
    void record_hospitalization(Agent<TSeq> & agent);

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
