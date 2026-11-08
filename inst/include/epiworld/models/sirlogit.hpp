#ifndef EPIWORLD_MODELS_SIRLOGIT_HPP 
#define EPIWORLD_MODELS_SIRLOGIT_HPP

#include "../model-bones.hpp"

/**
 * @brief A Susceptible-Infected-Removed (SIR) model
 * 
 * ![Model Diagram](../assets/img/sirlogit.png)
 * 
 * @ingroup special_models
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
 
*/
class ModelSIRLogit : public Model
{
public:
    static constexpr int SUSCEPTIBLE = 0;
    static constexpr int INFECTED    = 1;
    static constexpr int RECOVERED   = 2;

    ModelSIRLogit() = delete;

    /**
      * @param vname Name of the virus.
 * @param n Number of agents.
      * @param coefs_infect Double ptr. Infection coefficients.
      * @param coefs_recover Double ptr. Recovery coefficients.
      * @param ncoef_infect Unsigned int. Number of infection coefficients.
      * @param ncoef_recover Unsigned int. Number of recovery coefficients.
      * @param coef_infect_cols Vector<unsigned int>. Ids of infection vars.
      * @param coef_recover_cols Vector<unsigned int>. Ids of recover vars.
    */
    ModelSIRLogit(
        const std::string & vname,
        epiworld_fast_uint n,
        double * data,
        size_t ncols,
        std::vector< double > coefs_infect,
        std::vector< double > coefs_recover,
        std::vector< size_t > coef_infect_cols,
        std::vector< size_t > coef_recover_cols,
        epiworld_double transmission_prob,
        epiworld_double recovery_rate,
        epiworld_double prevalence
    );

    std::unique_ptr< Model > clone_ptr() override;

    void reset() override;
    
    std::vector< double > coefs_infect;
    std::vector< double > coefs_recover;
    std::vector< size_t > coef_infect_cols;
    std::vector< size_t > coef_recover_cols;

};

inline std::unique_ptr<Model> ModelSIRLogit::clone_ptr()
{
    
    return std::make_unique<ModelSIRLogit>(*this);

}

inline void ModelSIRLogit::reset()
{

    /* Checking specified columns in the model */
    for (const auto & c : coef_infect_cols)
    {
        if (c >= Model::agents_data_ncols)
            throw std::range_error("Columns specified in coef_infect_cols out of range.");
    }

    for (const auto & c : coef_recover_cols)
    {
        if (c >= Model::agents_data_ncols)
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
    
    Model::reset();

    return;

}

/**
 * @brief A Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param prob_transmission Probability of transmission
 * @param prob_recovery Probability of recovery
 */
inline ModelSIRLogit::ModelSIRLogit(
    const std::string & vname,
    epiworld_fast_uint n,
    double * data,
    size_t ncols,
    std::vector< double > coefs_infect,
    std::vector< double > coefs_recover,
    std::vector< size_t > coef_infect_cols,
    std::vector< size_t > coef_recover_cols,
    epiworld_double transmission_prob,
    epiworld_double recovery_rate,
    epiworld_double prevalence
    )
{

    this->set_population_size(n);

    if (coef_infect_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_infect_cols.");

    if (coef_recover_cols.size() == 0u)
        throw std::logic_error("No columns specified for coef_recover_cols.");

    // Saving the variables
    this->set_agents_data(
        data, ncols
    );

    this->coefs_infect = coefs_infect;
    this->coefs_recover = coefs_recover;
    this->coef_infect_cols = coef_infect_cols;
    this->coef_recover_cols = coef_recover_cols;

    this->add_param(transmission_prob, "Transmission probability");
    this->add_param(recovery_rate, "Recovery rate");

    // Infection: each infectious neighbor is an exposure with probability
    // plogis(baseline + p * coef_exposure), where the baseline comes from the
    // agent's features and p is the usual per-contact probability
    const Value p_contact("Transmission probability");
    InfectionSampler sampler = [p_contact](
        Agent & p, Model & m
        ) -> Virus *
        {

            ModelSIRLogit * _m = model_cast<ModelSIRLogit>(&m);

            const double coef_exposure = _m->coefs_infect[0u];

            double baseline = 0.0;
            for (size_t k = 0u; k < _m->coef_infect_cols.size(); ++k)
                baseline += p(k, m) * _m->coefs_infect[k + 1u];

            size_t nviruses_tmp = 0u;

            // Every network layer, in the order the layers were added
            for (auto * net : m.get_network_layers())
            for (auto * neighbor : net->neighbors_view(static_cast< size_t >(p.get_id())))
            {

                if (neighbor->get_virus() == nullptr)
                    continue;

                auto & v = neighbor->get_virus();

                m.reserve_tmp(nviruses_tmp + 1u);

                const double x = baseline +
                    (1.0 - p.get_susceptibility_reduction(v, m)) *
                    p_contact.eval(p, v.get(), m) *
                    (1.0 - neighbor->get_transmission_reduction(v, m)) *
                    coef_exposure;

                m.array_double_tmp[nviruses_tmp] = 1.0 / (1.0 + std::exp(-x));
                m.array_virus_tmp[nviruses_tmp++] = &(*v);

            }

            if (nviruses_tmp == 0u)
                return nullptr;

            const int which = roulette(nviruses_tmp, &m);
            return (which < 0) ? nullptr : m.array_virus_tmp[which];

        };

    // Recovery: plogis of the agent's features
    ValueFun p_recover = [](Agent & p, Model & m) -> epiworld_double
    {

        ModelSIRLogit * _m = model_cast<ModelSIRLogit>(&m);

        double prob = 0.0;
        for (size_t i = 0u; i < _m->coefs_recover.size(); ++i)
            prob += p(i, m) * _m->coefs_recover[i];

        return 1.0 / (1.0 + std::exp(-prob));

    };

    Virus virus(vname);
    virus.add_stage("Susceptible");
    virus.add_stage("Infected", stage::infectious);
    virus.add_stage("Recovered");
    virus.infection(SUSCEPTIBLE, INFECTED, "Transmission probability").sampler(sampler);
    virus.rate(INFECTED, RECOVERED, p_recover);

    this->add_virus(virus, assign::share(prevalence));

    this->set_name("Susceptible-Infected-Removed (SIR) (logit)");

}


#endif
