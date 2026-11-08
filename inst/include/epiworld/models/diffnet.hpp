#ifndef EPIWORLD_DIFFNET_H 
#define EPIWORLD_DIFFNET_H

#include <string>
#include <vector>
#include "../config.hpp"
#include "../model-bones.hpp"

/**
 * @brief A Network Diffusion Model
 * @ingroup special_models
 * 
 * @param model A Model object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param n Number of agents.
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery rate of the immune system
 * 
 */
class ModelDiffNet : public Model
{
private:
public:

    // Statuses
    static constexpr int NONADOPTER = 0;
    static constexpr int ADOPTER    = 1;

    ModelDiffNet() = delete;

    ModelDiffNet(
        const std::string & innovation_name,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double prob_adopt,
        bool normalize_exposure = true,
        double * agents_data = nullptr,
        size_t data_ncols = 0u,
        std::vector< size_t > data_cols = {},
        std::vector< double > params = {}
    );

    bool normalize_exposure = true;
    std::vector< size_t > data_cols;
    std::vector< double > params;
};

inline ModelDiffNet::ModelDiffNet(
    const std::string & innovation_name,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double prob_adopt,
    bool normalize_exposure,
    double * agents_data,
    size_t data_ncols,
    std::vector< size_t > data_cols,
    std::vector< double > params
    )
{

    this->set_population_size(n);

    // Adding additional parameters
    this->normalize_exposure = normalize_exposure;
    this->data_cols = data_cols;
    this->params = params;

    std::string parname = std::string("Prob. Adopting ") + innovation_name;
    this->add_param(prob_adopt, parname);

    // Adoption: for each innovation, a logit of the exposure (the share of
    // neighbors who adopted it, optionally normalized) plus the agent's
    // features and the baseline probability; one draw picks one, or none
    InfectionSampler sampler = [](
        Agent & agent, Model & m
    ) -> Virus * {

        const size_t nviruses = m.get_n_viruses();
        std::vector< Virus* > innovations(nviruses, nullptr);
        std::vector< double > exposure(nviruses, 0.0);

        ModelDiffNet * diffmodel = model_cast<ModelDiffNet>(&m);

        size_t n_neighbors = 0u;
        // Every network layer, in the order the layers were added
        for (auto * net : m.get_network_layers())
        for (auto * neighbor : net->neighbors_view(static_cast< size_t >(agent.get_id())))
        {

            ++n_neighbors;

            auto & v = neighbor->get_virus();
            if (v == nullptr)
                continue;

            const size_t vid = static_cast< size_t >(v->get_id());
            if (innovations[vid] == nullptr)
                innovations[vid] = &(*v);

            exposure[vid] +=
                (1.0 - agent.get_susceptibility_reduction(v, m)) *
                (1.0 - agent.get_transmission_reduction(v, m));

        }

        for (size_t i = 0u; i < nviruses; ++i)
        {

            if (diffmodel->normalize_exposure && (n_neighbors > 0u))
                exposure[i] /= static_cast< double >(n_neighbors);

            for (auto & j: diffmodel->data_cols)
                exposure[i] += agent(j, m) * diffmodel->params.at(j);

            // Baseline probability of adoption
            const double p = m.par(
                std::string("Prob. Adopting ") + m.get_viruses()[i]->get_name()
            );
            exposure[i] += std::log(p) - std::log(1.0 - p);

            exposure[i] = 1.0 / (1.0 + std::exp(-exposure[i]));

        }

        const int which = roulette(exposure, &m);
        if (which < 0)
            return nullptr;

        // An innovation no neighbor has is adopted from the model's own
        Virus * v = innovations[static_cast< size_t >(which)];
        return (v != nullptr) ? v : m.get_viruses()[static_cast< size_t >(which)].get();

    };

    this->set_agents_data(agents_data, data_ncols);

    Virus virus(innovation_name);
    virus.add_stage("Non adopters");
    virus.add_stage("Adopters", stage::infectious);
    virus.infection(NONADOPTER, ADOPTER, parname).sampler(sampler);

    this->add_virus(virus, assign::share(prevalence));

    this->set_name(
        std::string("Diffusion of Innovations - ") + innovation_name);

    return;
   
}

#endif
