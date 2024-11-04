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
 * @param params_now Vector where to save the new parameters.
 * @param params_prev Vector of reference parameters.
 * @param m LFMCMC model.
 * @tparam TData 
 */
template<typename TData>
inline void proposal_fun_normal(
    std::vector< epiworld_double >& params_now,
    const std::vector< epiworld_double >& params_prev,
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
 * @param params_now Where to write the new parameters
 * @param params_prev Reference parameters
 * @tparam TData 
 * @param m LFMCMC model.
 */
template<typename TData>
inline void proposal_fun_unif(
    std::vector< epiworld_double >& params_now,
    const std::vector< epiworld_double >& params_prev,
    LFMCMC<TData>* m
);

/**
 * @brief Uses the uniform kernel with euclidean distance
 * 
 * @param stats_now Vector of current statistics based on 
 * simulated data.
 * @param stats_obs Vector of observed statistics
 * @param epsilon Epsilon parameter
 * @param m LFMCMC model.
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_uniform(
    const std::vector< epiworld_double >& stats_now,
    const std::vector< epiworld_double >& stats_obs,
    epiworld_double epsilon,
    LFMCMC<TData>* m
);

/**
 * @brief Gaussian kernel
 * 
 * @tparam TData 
 * @param epsilon 
 * @param m 
 * @return epiworld_double 
 */
template<typename TData>
inline epiworld_double kernel_fun_gaussian(
    const std::vector< epiworld_double >& stats_now,
    const std::vector< epiworld_double >& stats_obs,
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
    std::shared_ptr< std::mt19937 > engine = nullptr;
    
    std::shared_ptr< std::uniform_real_distribution<> > runifd =
        std::make_shared< std::uniform_real_distribution<> >(0.0, 1.0);

    std::shared_ptr< std::normal_distribution<> > rnormd =
        std::make_shared< std::normal_distribution<> >(0.0);

    std::shared_ptr< std::gamma_distribution<> > rgammad = 
        std::make_shared< std::gamma_distribution<> >();

    // Process data
    TData observed_data;
    
    // Information about the size of the problem
    size_t n_samples;
    size_t n_statistics;
    size_t n_parameters;

    epiworld_double epsilon;

    std::vector< epiworld_double > params_now;
    std::vector< epiworld_double > params_prev;
    std::vector< epiworld_double > params_init;

    std::vector< epiworld_double > observed_stats; ///< Observed statistics

    std::vector< epiworld_double > sampled_params;     ///< Sampled Parameters
    std::vector< epiworld_double > sampled_stats;      ///< Sampled statistics
    std::vector< epiworld_double > sampled_stats_prob; ///< Sampled statistics
    std::vector< bool >            sampled_accepted;   ///< Indicator of accepted statistics

    std::vector< epiworld_double > accepted_params;      ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_stats;       ///< Posterior distribution (accepted samples)
    std::vector< epiworld_double > accepted_params_prob; ///< Posterior probability

    std::vector< epiworld_double > drawn_prob;     ///< Drawn probabilities (runif())
    std::vector< TData > * sampled_data = nullptr;

    // Functions
    LFMCMCSimFun<TData> simulation_fun;
    LFMCMCSummaryFun<TData> summary_fun;
    LFMCMCProposalFun<TData> proposal_fun = proposal_fun_normal<TData>;
    LFMCMCKernelFun<TData> kernel_fun     = kernel_fun_uniform<TData>;

    // Misc
    std::vector< std::string > names_parameters;
    std::vector< std::string > names_statistics;

    std::chrono::time_point<std::chrono::steady_clock> time_start;
    std::chrono::time_point<std::chrono::steady_clock> time_end;

    // std::chrono::milliseconds
    std::chrono::duration<epiworld_double,std::micro> time_elapsed = 
        std::chrono::duration<epiworld_double,std::micro>::zero();

    inline void get_elapsed(
        std::string unit,
        epiworld_double * last_elapsed,
        std::string * unit_abbr,
        bool print
    );

    void chrono_start();
    void chrono_end();
    
public:

    void run(
        std::vector< epiworld_double > param_init,
        size_t n_samples_,
        epiworld_double epsilon_,
        int seed = -1
        );

    LFMCMC() {};
    LFMCMC(const TData & observed_data_) : observed_data(observed_data_) {};
    ~LFMCMC() {};

    void set_observed_data(const TData & observed_data_) {observed_data = observed_data_;};
    void set_proposal_fun(LFMCMCProposalFun<TData> fun);
    void set_simulation_fun(LFMCMCSimFun<TData> fun);
    void set_summary_fun(LFMCMCSummaryFun<TData> fun);
    void set_kernel_fun(LFMCMCKernelFun<TData> fun);
    
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
    size_t get_n_samples() const {return n_samples;};
    size_t get_n_statistics() const {return n_statistics;};
    size_t get_n_parameters() const {return n_parameters;};
    epiworld_double get_epsilon() const {return epsilon;};

    const std::vector< epiworld_double > & get_params_now() {return params_now;};
    const std::vector< epiworld_double > & get_params_prev() {return params_prev;};
    const std::vector< epiworld_double > & get_params_init() {return params_init;};
    const std::vector< epiworld_double > & get_statistics_obs() {return observed_stats;};
    const std::vector< epiworld_double > & get_statistics_hist() {return sampled_stats;};
    const std::vector< bool >            & get_statistics_accepted() {return sampled_accepted;};
    const std::vector< epiworld_double > & get_posterior_lf_prob() {return accepted_params_prob;};
    const std::vector< epiworld_double > & get_drawn_prob() {return drawn_prob;};
    std::vector< TData > * get_sampled_data() {return sampled_data;};

    void set_par_names(std::vector< std::string > names);
    void set_stats_names(std::vector< std::string > names);

    std::vector< epiworld_double > get_params_mean();
    std::vector< epiworld_double > get_stats_mean();

    void print() ;

};

#endif