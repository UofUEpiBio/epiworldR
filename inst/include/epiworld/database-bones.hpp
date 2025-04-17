#ifndef EPIWORLD_DATABASE_BONES_HPP
#define EPIWORLD_DATABASE_BONES_HPP

template<typename TSeq>
class Model;

template<typename TSeq>
class Virus;

template<typename TSeq>
class UserData;

template<typename TSeq>
inline void default_add_virus(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_add_tool(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_virus(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_rm_tool(Event<TSeq> & a, Model<TSeq> * m);

template<typename TSeq>
inline void default_change_state(Event<TSeq> & a, Model<TSeq> * m);

/**
 * @brief Statistical data about the process
 * 
 * @tparam TSeq 
 */
template<typename TSeq>
class DataBase {
    friend class Model<TSeq>;
    friend void default_add_virus<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_add_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_virus<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_rm_tool<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
    friend void default_change_state<TSeq>(Event<TSeq> & a, Model<TSeq> * m);
private:
    Model<TSeq> * model;

    // Variants information 
    MapVec_type<int,int> virus_id; ///< The squence is the key
    std::vector< std::string > virus_name;
    std::vector< TSeq> virus_sequence;
    std::vector< int > virus_origin_date;
    std::vector< int > virus_parent_id;

    MapVec_type<int,int> tool_id; ///< The squence is the key
    std::vector< std::string > tool_name;
    std::vector< TSeq> tool_sequence;
    std::vector< int > tool_origin_date;

    std::function<std::vector<int>(const TSeq&)> seq_hasher = default_seq_hasher<TSeq>;
    std::function<std::string(const TSeq &)> seq_writer = default_seq_writer<TSeq>;

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_virus;

    // {Variant 1: {state 1, state 2, etc.}, Variant 2: {...}, ...}
    std::vector< std::vector<int> > today_tool;

    // {Susceptible, Infected, etc.}
    std::vector< int > today_total;

    // Totals
    int today_total_nviruses_active = 0;
    
    int sampling_freq = 1;

    // Variants history
    std::vector< int > hist_virus_date;
    std::vector< int > hist_virus_id;
    std::vector< epiworld_fast_uint > hist_virus_state;
    std::vector< int > hist_virus_counts;

    // Tools history
    std::vector< int > hist_tool_date;
    std::vector< int > hist_tool_id;
    std::vector< epiworld_fast_uint > hist_tool_state;
    std::vector< int > hist_tool_counts;

    // Overall hist
    std::vector< int > hist_total_date;
    std::vector< int > hist_total_nviruses_active;
    std::vector< epiworld_fast_uint > hist_total_state;
    std::vector< int > hist_total_counts;
    std::vector< int > hist_transition_matrix;

    // Transmission network
    std::vector< int > transmission_date;                 ///< Date of the transmission event
    std::vector< int > transmission_source;               ///< Id of the source
    std::vector< int > transmission_target;               ///< Id of the target
    std::vector< int > transmission_virus;              ///< Id of the variant
    std::vector< int > transmission_source_exposure_date; ///< Date when the source acquired the variant

    std::vector< int > transition_matrix;

    UserData<TSeq> user_data;

    void update_state(
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state,
        bool undo = false
    );

    void update_virus(
        epiworld_fast_uint virus_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
    );

    void update_tool(
        epiworld_fast_uint tool_id,
        epiworld_fast_uint prev_state,
        epiworld_fast_uint new_state
    );

    void record_transition(epiworld_fast_uint from, epiworld_fast_uint to, bool undo);


public:

    #ifdef EPI_DEBUG
    int n_transmissions_potential = 0;
    int n_transmissions_today     = 0;
    #endif

    DataBase() = delete;
    DataBase(Model<TSeq> & m) : model(&m), user_data(m) {};
    DataBase(const DataBase<TSeq> & db);
    // DataBase<TSeq> & operator=(const DataBase<TSeq> & m);

    /**
     * @brief Registering a new variant
     * 
     * @param v Pointer to the new virus.
     * Since viruses are originated in the agent, the numbers simply move around.
     * From the parent virus to the new virus. And the total number of infected
     * does not change.
     */
    void record_virus(Virus<TSeq> & v); 
    void record_tool(Tool<TSeq> & t); 
    void set_seq_hasher(std::function<std::vector<int>(TSeq)> fun);
    void reset();
    Model<TSeq> * get_model();
    void record();

    const std::vector< TSeq > & get_sequence() const;
    const std::vector< int > & get_nexposed() const;
    size_t size() const;

    /**
     * @name Get recorded information from the model
     * 
     * @param what std::string, The state, e.g., 0, 1, 2, ...
     * @return In `get_today_total`, the current counts of `what`.
     * @return In `get_today_virus`, the current counts of `what` for
     * each virus.
     * @return In `get_hist_total`, the time series of `what`
     * @return In `get_hist_virus`, the time series of what for each virus.
     * @return In `get_hist_total_date` and `get_hist_virus_date` the
     * corresponding date
     */
    ///@{
    int get_today_total(std::string what) const;
    int get_today_total(epiworld_fast_uint what) const;
    void get_today_total(
        std::vector< std::string > * state = nullptr,
        std::vector< int > * counts = nullptr
    ) const;

    void get_today_virus(
        std::vector< std::string > & state,
        std::vector< int > & id,
        std::vector< int > & counts
    ) const;

    void get_today_transition_matrix(
        std::vector< int > & counts
    ) const;

    void get_hist_total(
        std::vector< int > * date,
        std::vector< std::string > * state,
        std::vector< int > * counts
    ) const;

    void get_hist_virus(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & state,
        std::vector< int > & counts
    ) const;

    void get_hist_tool(
        std::vector< int > & date,
        std::vector< int > & id,
        std::vector< std::string > & state,
        std::vector< int > & counts
    ) const;

    void get_hist_transition_matrix(
        std::vector< std::string > & state_from,
        std::vector< std::string > & state_to,
        std::vector< int > & date,
        std::vector< int > & counts,
        bool skip_zeros
    ) const;
    ///@}

    /**
     * @brief Get the transmissions object
     * 
     * @param date 
     * @param source 
     * @param target 
     * @param virus 
     * @param source_exposure_date 
     */
    ///@{
    void get_transmissions(
        std::vector<int> & date,
        std::vector<int> & source,
        std::vector<int> & target,
        std::vector<int> & virus,
        std::vector<int> & source_exposure_date
    ) const;

    void get_transmissions(
        int * date,
        int * source,
        int * target,
        int * virus,
        int * source_exposure_date
    ) const;
    ///@}

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
    
    /***
     * @brief Record a transmission event
     * @param i,j Integers. Id of the source and target agents.
     * @param virus Integer. Id of the virus.
     * @param i_expo_date Integer. Date when the source agent was infected.
     * @details
     * If i is -1, then it means that the agent was assigned a virus at the
     * beginning of the simulation.
     */
    void record_transmission(int i, int j, int virus, int i_expo_date);

    size_t get_n_viruses() const; ///< Get the number of viruses
    size_t get_n_tools() const; ///< Get the number of tools
    
    void set_user_data(std::vector< std::string > names);
    void add_user_data(std::vector< epiworld_double > x);
    void add_user_data(epiworld_fast_uint j, epiworld_double x);
    UserData<TSeq> & get_user_data();


    /**
     * @brief Computes the reproductive number of each case
     * 
     * @details By definition, whereas it computes R0 (basic reproductive number)
     * or Rt/R (the effective reproductive number) will depend on whether the
     * virus is allowed to circulate naïvely or not, respectively.
     * 
     * @param fn File where to write out the reproductive number.
     * @details
     * In the case of `MapVec_type<int,int>`, the key is a vector of 3 integers:
     * - Virus id
     * - Source id
     * - Date when the source was infected
     */
    ///@{
    MapVec_type<int,int> reproductive_number() const;

    void reproductive_number(
        std::string fn
        ) const;
    ///@}

    /**
     * @brief Calculates the transition probabilities
     * @param print Print the transition matrix.
     * @param normalize Normalize the transition matrix. Otherwise, 
     * it returns raw counts.
     * @details
     * The transition matrix is the matrix of the counts of transitions
     * from one state to another. So the ij-th element of the matrix is
     * the number of transitions from state i to state j (when not normalized),
     * or the probability of transitioning from state i to state j 
     * (when normalized).
     * @return std::vector< epiworld_double > 
     */
    std::vector< epiworld_double > transition_probability(
        bool print = true,
        bool normalize = true
    ) const;

    bool operator==(const DataBase<TSeq> & other) const;
    bool operator!=(const DataBase<TSeq> & other) const {return !operator==(other);};

    /**
     * Calculates the generating time
     * @param agent_id,virus_id,time,gentime vectors where to save the values 
     * 
     * @details
     * The generation time is the time between the infection of the source and 
     * the infection of the target.
    */
   ///@{
    void generation_time(
        std::vector< int > & agent_id,
        std::vector< int > & virus_id,
        std::vector< int > & time,
        std::vector< int > & gentime
    ) const; ///< Get the generation time

    void generation_time(
        std::string fn
    ) const; ///< Write the generation time to a file
    ///@}

};


#endif