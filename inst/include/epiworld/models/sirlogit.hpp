// #include "../epiworld.hpp"

#ifndef EPIWORLD_MODELS_SIRLOGIT_HPP 
#define EPIWORLD_MODELS_SIRLOGIT_HPP


/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @details
 * In this model, infection and recoveru probabilities are computed
 * using a logit model. Particularly, the probability of infection
 * is computed as:
 * 
 * \f[
 * \frac{1}{1 + \exp\left(-\left(\beta_0 E_i + \sum_{i=1}^{n} \beta_i x_i\right)\right)}
 * \f]
 * 
 * where \f$\beta_0\f$ is the exposure coefficient and \f$E_i\f$ is the exposure
 * number, \f$\beta_i\f$ are the
 * coefficients for the features \f$x_i\f$ of the agents, and \f$n\f$ is the
 * number of features. The probability of recovery is computed as:
 * 
 * \f[
 * \frac{1}{1 + \exp\left(-\left(\sum_{i=1}^{n} \beta_i x_i\right)\right)}
 * \f]
 * 
 * where \f$\beta_i\f$ are the coefficients for the features \f$x_i\f$ of the agents,
 * and \f$n\f$ is the number of features.
 * 
 * @param TSeq Type of the sequence (e.g. std::vector, std::deque)
 
*/
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRLogit : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;

public:

    ModelSIRLogit() {};

    /**
      * @param vname Name of the virus.
      * @param coefs_infect Double ptr. Infection coefficients.
      * @param coefs_recover Double ptr. Recovery coefficients.
      * @param ncoef_infect Unsigned int. Number of infection coefficients.
      * @param ncoef_recover Unsigned int. Number of recovery coefficients.
      * @param coef_infect_cols Vector<unsigned int>. Ids of infection vars.
      * @param coef_recover_cols Vector<unsigned int>. Ids of recover vars.
    */
    ModelSIRLogit(
        ModelSIRLogit<TSeq> & model,
        const std::string & vname,
        double * data,
        size_t ncols,
        std::vector< double > coefs_infect,
        std::vector< double > coefs_recover,
        std::vector< size_t > coef_infect_cols,
        std::vector< size_t > coef_recover_cols,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate,
        epiworld_double prevalence
    );

    ModelSIRLogit(
        const std::string & vname,
        double * data,
        size_t ncols,
        std::vector< double > coefs_infect,
        std::vector< double > coefs_recover,
        std::vector< size_t > coef_infect_cols,
        std::vector< size_t > coef_recover_cols,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate,
        epiworld_double prevalence
    );

    ModelSIRLogit<TSeq> & run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    Model<TSeq> * clone_ptr();

    void reset();
    
    std::vector< double > coefs_infect;
    std::vector< double > coefs_recover;
    std::vector< size_t > coef_infect_cols;
    std::vector< size_t > coef_recover_cols;

};



template<typename TSeq>
inline ModelSIRLogit<TSeq> & ModelSIRLogit<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    Model<TSeq>::run(ndays, seed);
    return *this;

}

template<typename TSeq>
inline Model<TSeq> * ModelSIRLogit<TSeq>::clone_ptr()
{
    
    ModelSIRLogit<TSeq> * ptr = new ModelSIRLogit<TSeq>(
        *dynamic_cast<const ModelSIRLogit<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

template<typename TSeq>
inline void ModelSIRLogit<TSeq>::reset()
{

    /* Checking specified columns in the model */
    for (const auto & c : coef_infect_cols)
    {
        if (c >= Model<TSeq>::agents_data_ncols)
            throw std::range_error("Columns specified in coef_infect_cols out of range.");
    }

    for (const auto & c : coef_recover_cols)
    {
        if (c >= Model<TSeq>::agents_data_ncols)
            throw std::range_error("Columns specified in coef_recover_cols out of range.");
    }

    /* Checking attributes */ 
    if (coefs_infect.size() != (coef_infect_cols.size() + 1u))
        throw std::logic_error(
            "The number of coefficients (infection) doesn't match the number of features. It must be as many features of the agents plus 1 (exposure.)"
            );

    if (coefs_recover.size() != coef_recover_cols.size())
        throw std::logic_error(
            "The number of coefficients (recovery) doesn't match the number of features. It must be as many features of the agents."
            );
    
    Model<TSeq>::reset();

    return;

}

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param prob_transmission Probability of transmission
 * @param prob_recovery Probability of recovery
 */
template<typename TSeq>
inline ModelSIRLogit<TSeq>::ModelSIRLogit(
    ModelSIRLogit<TSeq> & model,
    const std::string & vname,
    double * data,
    size_t ncols,
    std::vector< double > coefs_infect,
    std::vector< double > coefs_recover,
    std::vector< size_t > coef_infect_cols,
    std::vector< size_t > coef_recover_cols,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double prevalence
    )
{

    if (coef_infect_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_infect_cols.");

    if (coef_recover_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_recover_cols.");

    // Saving the variables
    model.set_agents_data(
        data, ncols
    );

    model.coefs_infect = coefs_infect;
    model.coefs_recover = coefs_recover;
    model.coef_infect_cols = coef_infect_cols;
    model.coef_recover_cols = coef_recover_cols;

    epiworld::UpdateFun<TSeq> update_susceptible = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRLogit<TSeq> * _m = dynamic_cast<ModelSIRLogit<TSeq>*>(m);

            // Exposure coefficient
            const double coef_exposure = _m->coefs_infect[0u];

            // This computes the prob of getting any neighbor variant
            size_t nviruses_tmp = 0u;

            double baseline = 0.0;
            for (size_t k = 0u; k < _m->coef_infect_cols.size(); ++k)
                baseline += p->operator[](k) * _m->coefs_infect[k + 1u];

            for (auto & neighbor: p->get_neighbors()) 
            {
                
                if (neighbor->get_virus() == nullptr)
                    continue;

                auto & v = neighbor->get_virus();

                #ifdef EPI_DEBUG
                if (nviruses_tmp >= m->array_virus_tmp.size())
                    throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                #endif
                    
                /* And it is a function of susceptibility_reduction as well */ 
                m->array_double_tmp[nviruses_tmp] =
                    baseline +
                    (1.0 - p->get_susceptibility_reduction(v, m)) * 
                    v->get_prob_infecting(m) * 
                    (1.0 - neighbor->get_transmission_reduction(v, m))  *
                    coef_exposure
                    ; 

                // Applying the plogis function
                m->array_double_tmp[nviruses_tmp] = 1.0/
                    (1.0 + std::exp(-m->array_double_tmp[nviruses_tmp]));
            
                m->array_virus_tmp[nviruses_tmp++] = &(*v);

            }

            // No virus to compute
            if (nviruses_tmp == 0u)
                return;

            // Running the roulette
            int which = roulette(nviruses_tmp, m);

            if (which < 0)
                return;

            p->set_virus(*m->array_virus_tmp[which], m);

            return;

        };

    epiworld::UpdateFun<TSeq> update_infected = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Getting the right type
            ModelSIRLogit<TSeq> * _m = dynamic_cast<ModelSIRLogit<TSeq>*>(m);

            // Computing recovery probability once
            double prob    = 0.0;
            #if defined(__OPENMP) || defined(_OPENMP)
            #pragma omp simd reduction(+:prob)
            #endif
            for (size_t i = 0u; i < _m->coefs_recover.size(); ++i)
                prob += p->operator[](i) * _m->coefs_recover[i];

            // Computing logis
            prob = 1.0/(1.0 + std::exp(-prob));

            if (prob > m->runif())
                p->rm_virus(m);
            
            return;

        };

    // state
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Setting up parameters
    // model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(recovery_rate, "Recovery rate");
    // model.add_param(prob_reinfection, "Prob. Reinfection");
    
    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(
        ModelSIRLogit<TSeq>::INFECTED,
        ModelSIRLogit<TSeq>::RECOVERED,
        ModelSIRLogit<TSeq>::RECOVERED
        );

    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_prob_recovery(&model("Recovery rate"));

    // virus.set_prob

    model.add_virus(virus);

    model.set_name("Susceptible-Infected-Removed (SIR) (logit)");

    return;

}

template<typename TSeq>
inline ModelSIRLogit<TSeq>::ModelSIRLogit(
    const std::string & vname,
    double * data,
    size_t ncols,
    std::vector< double > coefs_infect,
    std::vector< double > coefs_recover,
    std::vector< size_t > coef_infect_cols,
    std::vector< size_t > coef_recover_cols,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double prevalence
    )
{

    ModelSIRLogit(
        *this,
        vname,
        data,
        ncols,
        coefs_infect,
        coefs_recover,
        coef_infect_cols,
        coef_recover_cols,
        transmission_rate,
        recovery_rate,
        prevalence
    );

    return;

}


#endif
