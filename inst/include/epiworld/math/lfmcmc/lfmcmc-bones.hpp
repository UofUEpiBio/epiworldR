#ifndef EPIWORLD_LFMCMC_BONES_HPP
#define EPIWORLD_LFMCMC_BONES_HPP

#ifndef epiworld_double
    #define epiworld_double float
#endif


template<typename TData>
class LFMCMC;

template<typename TData>
using LFMCMCSimFun = std::function<TData(const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCSummaryFun = std::function<void(std::vector< epiworld_double >&,const TData&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCProposalFun = std::function<void(std::vector< epiworld_double >&,const std::vector< epiworld_double >&,LFMCMC<TData>*)>;

template<typename TData>
using LFMCMCKernelFun = std::function<epiworld_double(const std::vector< epiworld_double >&,const std::vector< epiworld_double >&,epiworld_double,LFMCMC<TData>*)>;

/**
 * @brief Proposal function
 * @param new_params Vector where to save the new parameters.
 * @param old_params Vector of reference parameters.
 * @param m LFMCMC model.
 * @tparam TData 
 */
template<typename TData>
inline void proposal_fun_normal(
    std::vector< epiworld_double >& new_params,
    const std::vector< epiworld_double >& old_params,
    LFMCMC<TData>* m
);

/**
 * @brief Factory for a reflective normal kernel
 * 
 * @details Reflective kernel corrects proposals by forcing them to be
 * within prespecified boundaries. 
 * 
 * @tparam TData 
 * @param scale Scale of the normal kernel
 * @param lb Lower bound (applies the same to all parameters)
 * @param ub Upper bound (applies the same to all parameters)
 * @return LFMCMCProposalFun<TData> 
 */
template<typename TData>
inline LFMCMCProposalFun<TData> make_proposal_norm_reflective(
    epiworld_double scale,
    epiworld_double lb = std::numeric_limits<epiworld_double>::min(),
    epiworld_double ub = std::numeric_limits<epiworld_double>::max()
);

/**
 * @brief Uniform proposal kernel
 * 
 * Proposals are made within a radious 1 of the current
 * state of the parameters.
 * 
 * @param new_params Where to write the new parameters
 * @param old_params Reference parameters
 * @tparam TData 
 * @param m LFMCMC model.
 */
template<typename TData>
inline void proposal_fun_unif(
    std::vector< epiworld_double >& new_params,
    const std::vector< epiworld_double >& old_params,
    LFMCMC<TData>* m
);

/**
 * @brief Uses the uniform kernel with euclidean distance
 * 
 * @param simulated_stats Vector of statistics based on 
 * simulated data
 * @param observed_stats Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_uniform(
    const std::vector< epiworld_double >& simulated_stats,
    const std::vector< epiworld_double >& observed_stats,
    epiworld_double epsilon,
    LFMCMC<TData>* m
);

/**
 * @brief Gaussian kernel
 * 
 * @tparam TData 
 * @param simulated_stats Vector of statistics based on 
 * simulated data
 * @param observed_stats Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_gaussian(
    const std::vector< epiworld_double >& simulated_stats,
    const std::vector< epiworld_double >& observed_stats,
    epiworld_double epsilon,
    LFMCMC<TData>* m
);

/**
 * @brief Likelihood-Free Markov Chain Monte Carlo
 * 
 * @tparam TData Type of data that is generated
 */
template<typename TData>
class LFMCMC {
private:

    // Random number sampling
    std::shared_ptr< std::mt19937 > m_engine = nullptr;
    
    std::shared_ptr< std::uniform_real_distribution<> > runifd =
        std::make_shared< std::uniform_real_distribution<> >(0.0, 1.0);

    std::shared_ptr< std::normal_distribution<> > rnormd =
        std::make_shared< std::normal_distribution<> >(0.0);

    std::shared_ptr< std::gamma_distribution<> > rgammad = 
        std::make_shared< std::gamma_distribution<> >();

    // Process data
    TData m_observed_data;
    std::vector< TData > * m_simulated_data = nullptr;

    // Information about the size of the process
    size_t m_n_samples;
    size_t m_n_stats;
    size_t m_n_params;

    epiworld_double m_epsilon;

    std::vector< epiworld_double > m_initial_params;             ///< Initial parameters
    std::vector< epiworld_double > m_current_proposed_params;    ///< Proposed parameters for the next sample
    std::vector< epiworld_double > m_current_accepted_params;    ///< Most recently accepted parameters (current state of MCMC)
    std::vector< epiworld_double > m_current_proposed_stats;     ///< Statistics from simulation run with proposed parameters
    std::vector< epiworld_double > m_current_accepted_stats;     ///< Statistics from simulation run with most recently accepted params

    std::vector< epiworld_double > m_observed_stats;             ///< Observed statistics

    std::vector< epiworld_double > m_all_sample_params;          ///< Parameter samples
    std::vector< epiworld_double > m_all_sample_stats;           ///< Statistic samples
    std::vector< bool >            m_all_sample_acceptance;      ///< Indicator if sample was accepted
    std::vector< epiworld_double > m_all_sample_drawn_prob;      ///< Drawn probabilities (runif()) for each sample
    std::vector< epiworld_double > m_all_sample_kernel_scores;   ///< Kernel scores for each sample

    std::vector< epiworld_double > m_all_accepted_params;        ///< Posterior distribution of parameters from accepted samples
    std::vector< epiworld_double > m_all_accepted_stats;         ///< Posterior distribution of statistics from accepted samples
    std::vector< epiworld_double > m_all_accepted_kernel_scores; ///< Kernel scores for each accepted sample

    // Functions
    LFMCMCSimFun<TData> m_simulation_fun;
    LFMCMCSummaryFun<TData> m_summary_fun;
    LFMCMCProposalFun<TData> m_proposal_fun = proposal_fun_normal<TData>;
    LFMCMCKernelFun<TData> m_kernel_fun     = kernel_fun_uniform<TData>;

    // Misc
    std::vector< std::string > m_param_names;
    std::vector< std::string > m_stat_names;

    std::chrono::time_point<std::chrono::steady_clock> m_start_time;
    std::chrono::time_point<std::chrono::steady_clock> m_end_time;

    // Timing
    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> m_elapsed_time = 
        std::chrono::duration<epiworld_double,std::micro>::zero();

    inline void get_elapsed_time(
        std::string unit,
        epiworld_double * last_elapsed,
        std::string * unit_abbr,
        bool print
    ) const;

    void chrono_start();
    void chrono_end();

    // Progress
    bool verbose = true;
    Progress progress_bar;
    
public:

    void run(
        std::vector< epiworld_double > params_init_,
        size_t n_samples_,
        epiworld_double epsilon_,
        int seed = -1
        );

    LFMCMC() {};
    LFMCMC(const TData & observed_data_) : m_observed_data(observed_data_) {};
    ~LFMCMC() {};

    // Setting LFMCMC variables
    void set_observed_data(const TData & observed_data_) {m_observed_data = observed_data_;};
    
    void set_proposal_fun(LFMCMCProposalFun<TData> fun);
    void set_simulation_fun(LFMCMCSimFun<TData> fun);
    void set_summary_fun(LFMCMCSummaryFun<TData> fun);
    void set_kernel_fun(LFMCMCKernelFun<TData> fun);

    void set_params_names(std::vector< std::string > names);
    void set_stats_names(std::vector< std::string > names);
    
    /**
     * @name Random number generation
     * 
     * @param eng 
     */
    ///@{
    void set_rand_engine(std::shared_ptr< std::mt19937 > & eng);
    std::shared_ptr< std::mt19937 > & get_rand_endgine();
    void seed(epiworld_fast_uint s);
    void set_rand_gamma(epiworld_double alpha, epiworld_double beta);
    epiworld_double runif();
    epiworld_double rnorm();
    epiworld_double rgamma();
    epiworld_double runif(epiworld_double lb, epiworld_double ub);
    epiworld_double rnorm(epiworld_double mean, epiworld_double sd);
    epiworld_double rgamma(epiworld_double alpha, epiworld_double beta);
    ///@}

    // Accessing parameters of the function
    size_t get_n_samples() const {return m_n_samples;};
    size_t get_n_stats() const {return m_n_stats;};
    size_t get_n_params() const {return m_n_params;};
    epiworld_double get_epsilon() const {return m_epsilon;};

    const std::vector< epiworld_double > & get_initial_params() const {return m_initial_params;};
    const std::vector< epiworld_double > & get_current_proposed_params() const {return m_current_proposed_params;};
    const std::vector< epiworld_double > & get_current_accepted_params() const {return m_current_accepted_params;};
    const std::vector< epiworld_double > & get_current_proposed_stats() const {return m_current_proposed_stats;};
    const std::vector< epiworld_double > & get_current_accepted_stats() const {return m_current_accepted_stats;};

    const std::vector< epiworld_double > & get_observed_stats() const {return m_observed_stats;};

    const std::vector< epiworld_double > & get_all_sample_params() const {return m_all_sample_params;};
    const std::vector< epiworld_double > & get_all_sample_stats() const {return m_all_sample_stats;};
    const std::vector< bool >            & get_all_sample_acceptance() const {return m_all_sample_acceptance;};
    const std::vector< epiworld_double > & get_all_sample_drawn_prob() const {return m_all_sample_drawn_prob;};
    const std::vector< epiworld_double > & get_all_sample_kernel_scores() const {return m_all_sample_kernel_scores;};

    const std::vector< epiworld_double > & get_all_accepted_params() const {return m_all_accepted_params;};
    const std::vector< epiworld_double > & get_all_accepted_stats() const {return m_all_accepted_stats;};
    const std::vector< epiworld_double > & get_all_accepted_kernel_scores() const {return m_all_accepted_kernel_scores;};
    
    std::vector< TData > * get_simulated_data() const {return m_simulated_data;};

    std::vector< epiworld_double > get_mean_params();
    std::vector< epiworld_double > get_mean_stats();

    // Printing
    LFMCMC<TData> & verbose_off();
    LFMCMC<TData> & verbose_on();
    void print(size_t burnin = 0u) const;

};

#endif