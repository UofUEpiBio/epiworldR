#ifndef EPIWORLD_MODEL_BONES_HPP
#define EPIWORLD_MODEL_BONES_HPP

template<typename TSeq>
class Agent;

template<typename TSeq>
class AgentsSample;

template<typename TSeq>
class Virus;

template<typename TSeq>
class Viruses;

template<typename TSeq>
class Viruses_const;

template<typename TSeq>
class Tool;

class AdjList;

template<typename TSeq>
class DataBase;

template<typename TSeq>
class Queue;

template<typename TSeq>
struct Event;

template<typename TSeq>
class GlobalEvent;

template<typename TSeq>
inline epiworld_double susceptibility_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double transmission_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double recovery_enhancer_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );
template<typename TSeq>
inline epiworld_double death_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
    );

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
    bool generation = false
    );

// template<typename TSeq>
// class VirusPtr;

// template<typename TSeq>
// class ToolPtr;

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
    std::vector< Entity<TSeq> > entities_backup = {};

    std::shared_ptr< std::mt19937 > engine = std::make_shared< std::mt19937 >();

    std::uniform_real_distribution<> runifd      =
        std::uniform_real_distribution<> (0.0, 1.0);
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

    void dist_tools();
    void dist_virus();
    void dist_entities();

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed =
        std::chrono::duration<epiworld_double,std::micro>::zero();
    epiworld_fast_uint n_replicates = 0u;
    void chrono_start();
    void chrono_end();

    std::vector<GlobalEvent<TSeq>> globalevents;

    Queue<TSeq> queue;
    bool use_queuing   = true;

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
     * @param call_ Function the action will call
     * @param queue_ Change in the queue
     * @param idx_agent_ Location of agent in object.
     * @param idx_object_ Location of object in agent.
     */
    void events_add(
        Agent<TSeq> * agent_,
        VirusPtr<TSeq> virus_,
        ToolPtr<TSeq> tool_,
        Entity<TSeq> * entity_,
        epiworld_fast_int new_state_,
        epiworld_fast_int queue_,
        EventFun<TSeq> call_,
        int idx_agent_,
        int idx_object_
        );

    /**
     * @name Tool Mixers
     *
     * These functions combine the effects tools have to deliver
     * a single effect. For example, wearing a mask, been vaccinated,
     * and the immune system combine together to jointly reduce
     * the susceptibility for a given virus.
     *
     */
    MixerFun<TSeq> susceptibility_reduction_mixer = susceptibility_reduction_mixer_default<TSeq>;
    MixerFun<TSeq> transmission_reduction_mixer = transmission_reduction_mixer_default<TSeq>;
    MixerFun<TSeq> recovery_enhancer_mixer = recovery_enhancer_mixer_default<TSeq>;
    MixerFun<TSeq> death_reduction_mixer = death_reduction_mixer_default<TSeq>;

    /**
     * @brief Advanced usage: Makes a copy of data and returns it as undeleted pointer
     *
     * @param copy
     */
    virtual Model<TSeq> * clone_ptr();

public:


    std::vector<epiworld_double> array_double_tmp;
    std::vector<Virus<TSeq> * > array_virus_tmp;

    Model();
    Model(const Model<TSeq> & m);
    Model(Model<TSeq> & m);
    Model(Model<TSeq> && m);
    Model<TSeq> & operator=(const Model<TSeq> & m);

    virtual ~Model() {};

    /**
     * @name Set the backup object
     * @details `backup` can be used to restore the entire object
     * after a run. This can be useful if the user wishes to have
     * individuals start with the same network from the beginning.
     *
     */
    ///@{
    void set_backup();
    // void restore_backup();
    ///@}

    DataBase<TSeq> & get_db();
    const DataBase<TSeq> & get_db() const;
    epiworld_double & operator()(std::string pname);

    size_t size() const;

    /**
     * @name Random number generation
     *
     * @param eng Random number generator
     * @param s Seed
     */
    ///@{
    void set_rand_engine(std::shared_ptr< std::mt19937 > & eng);
    std::shared_ptr< std::mt19937 > & get_rand_endgine();
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
    epiworld_double runif();
    epiworld_double runif(epiworld_double a, epiworld_double b);
    epiworld_double rnorm();
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma();
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double rexp();
    epiworld_double rexp(epiworld_double lambda);
    epiworld_double rlognormal();
    epiworld_double rlognormal(epiworld_double mean, epiworld_double shape);
    int rbinom();
    int rbinom(int n, epiworld_double p);
    int rnbinom();
    int rnbinom(int n, epiworld_double p);
    int rgeom();
    int rgeom(epiworld_double p);
    int rpoiss();
    int rpoiss(epiworld_double lambda);
    ///@}

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
     * @param fn std::string Filename of the edgelist file.
     * @param skip int Number of lines to skip in `fn`.
     * @param directed bool Whether the graph is directed or not.
     * @param size Size of the network.
     * @param al AdjList to read into the model.
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

    bool is_directed() const;

    std::vector< Agent<TSeq> > & get_agents(); ///< Returns a reference to the vector of agents.

    Agent<TSeq> & get_agent(size_t i);

    std::vector< epiworld_fast_uint > get_agents_states() const; ///< Returns a vector with the states of the agents.

    std::vector< Viruses_const<TSeq> > get_agents_viruses() const; ///< Returns a const vector with the viruses of the agents.

    std::vector< Viruses<TSeq> > get_agents_viruses(); ///< Returns a vector with the viruses of the agents.

    std::vector< Entity<TSeq> > & get_entities();

    Entity<TSeq> & get_entity(size_t entity_id, int * entity_pos = nullptr);

    Model<TSeq> & agents_smallworld(
        epiworld_fast_uint n = 1000,
        epiworld_fast_uint k = 5,
        bool d = false,
        epiworld_double p = .01
        );
    void agents_empty_graph(epiworld_fast_uint n = 1000);
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
    void next();
    virtual Model<TSeq> & run(
        epiworld_fast_uint ndays,
        int seed = -1
    ); ///< Runs the simulation (after initialization)
    void run_multiple( ///< Multiple runs of the simulation
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
    epiworld_fast_uint get_n_replicates() const;
    size_t get_n_entities() const;
    void set_ndays(epiworld_fast_uint ndays);
    bool get_verbose() const;
    Model<TSeq> & verbose_off();
    Model<TSeq> & verbose_on();
    int today() const; ///< The current time of the model

    /**
     * @name Rewire the network preserving the degree sequence.
     *
     * @details This implementation assumes an undirected network,
     * thus if {(i,j), (k,l)} -> {(i,l), (k,j)}, the reciprocal
     * is also true, i.e., {(j,i), (l,k)} -> {(j,k), (l,i)}.
     *
     * @param proportion Proportion of ties to be rewired.
     *
     * @result A rewired version of the network.
     */
    ///@{
    void set_rewire_fun(std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> fun);
    void set_rewire_prop(epiworld_double prop);
    epiworld_double get_rewire_prop() const;
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
        std::string fn_generation_time
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
     * states included in the model.
     *
     * @param lab `std::string` Name of the state.
     *
     * @return `add_state*` returns nothing.
     * @return `get_state_*` returns a vector of pairs with the
     * states and their labels.
     */
    ///@{
    void add_state(std::string lab, UpdateFun<TSeq> fun = nullptr);
    const std::vector< std::string > & get_states() const;
    size_t get_n_states() const;
    const std::vector< UpdateFun<TSeq> > & get_state_fun() const;
    void print_state_codes() const;
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
    epiworld_double get_param(epiworld_fast_uint k);
    epiworld_double get_param(std::string pname);
    // void set_param(size_t k, epiworld_double val);
    void set_param(std::string pname, epiworld_double val);
    // epiworld_double par(epiworld_fast_uint k);
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
        GlobalEvent<TSeq> action
    );

    GlobalEvent<TSeq> & get_globalevent(std::string name); ///< Retrieve a global action by name
    GlobalEvent<TSeq> & get_globalevent(size_t i); ///< Retrieve a global action by index

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
     * @name Get the susceptibility reduction object
     *
     * @param v
     * @return epiworld_double
     */
    ///@{
    void set_susceptibility_reduction_mixer(MixerFun<TSeq> fun);
    void set_transmission_reduction_mixer(MixerFun<TSeq> fun);
    void set_recovery_enhancer_mixer(MixerFun<TSeq> fun);
    void set_death_reduction_mixer(MixerFun<TSeq> fun);
    ///@}

    const std::vector< VirusPtr<TSeq> > & get_viruses() const;
    const std::vector< ToolPtr<TSeq> > & get_tools() const;
    Virus<TSeq> & get_virus(size_t id);
    Tool<TSeq> & get_tool(size_t id);

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


};

#endif
