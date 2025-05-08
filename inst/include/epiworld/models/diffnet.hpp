#ifndef EPIWORLD_DIFFNET_H 
#define EPIWORLD_DIFFNET_H

/**
 * @brief Template for a Network Diffusion Model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param initial_prevalence epiworld_double Initial prevalence
 * @param initial_efficacy epiworld_double Initial susceptibility_reduction of the immune system
 * @param initial_recovery epiworld_double Initial recovery rate of the immune system
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelDiffNet : public epiworld::Model<TSeq>
{
private:
public:

    ModelDiffNet() {};

    ModelDiffNet(
        ModelDiffNet<TSeq> & model,
        const std::string & innovation_name,
        epiworld_double prevalence,
        epiworld_double prob_adopt,
        bool normalize_exposure = true,
        double * agents_data = nullptr,
        size_t data_ncols = 0u,
        std::vector< size_t > data_cols = {},
        std::vector< double > params = {}
    );

    ModelDiffNet(
        const std::string & innovation_name,
        epiworld_double prevalence,
        epiworld_double prob_adopt,
        bool normalize_exposure = true,
        double * agents_data = nullptr,
        size_t data_ncols = 0u,
        std::vector< size_t > data_cols = {},
        std::vector< double > params = {}
    );
    
    static const int NONADOPTER = 0;
    static const int ADOPTER    = 1;

    bool normalize_exposure = true;
    std::vector< size_t > data_cols;
    std::vector< double > params;
};

template<typename TSeq>
inline ModelDiffNet<TSeq>::ModelDiffNet(
    ModelDiffNet<TSeq> & model,
    const std::string & innovation_name,
    epiworld_double prevalence,
    epiworld_double prob_adopt,
    bool normalize_exposure,
    double * agents_data,
    size_t data_ncols,
    std::vector< size_t > data_cols,
    std::vector< double > params
    )
{

    // Adding additional parameters
    this->normalize_exposure = normalize_exposure;
    this->data_cols = data_cols;
    this->params = params;

    epiworld::UpdateFun<TSeq> update_non_adopters = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
    ) -> void {

        // Measuring exposure
        // If the neighbor is infected, then proceed
        size_t nviruses = m->get_n_viruses();
        std::vector< Virus<TSeq>* > innovations(nviruses, {});
        std::vector< bool > stored(nviruses, false);
        std::vector< double > exposure(nviruses, 0.0);

        ModelDiffNet<TSeq> * diffmodel = dynamic_cast<ModelDiffNet<TSeq>*>(m);

        Agent<TSeq> & agent = *p;

        // For each one of the possible innovations, we have to compute
        // the adoption probability, which is a function of exposure
        for (auto & neighbor: agent.get_neighbors())
        {

            if (neighbor->get_state() == ModelDiffNet<TSeq>::ADOPTER)
            {

                auto & v = neighbor->get_virus();
                
                if (v == nullptr)
                    continue;
    
                /* And it is a function of susceptibility_reduction as well */ 
                double p_i =
                    (1.0 - agent.get_susceptibility_reduction(v, m)) * 
                    (1.0 - agent.get_transmission_reduction(v, m)) 
                    ; 
            
                size_t vid = v->get_id();
                if (!stored[vid])
                {
                    stored[vid] = true;
                    innovations[vid] = &(*v);
                }
                exposure[vid] += p_i;


            }

        }

        // Computing probability of adoption
        for (size_t i = 0u; i < nviruses; ++i)
        {

            if (diffmodel->normalize_exposure)
                exposure.at(i) /= agent.get_n_neighbors();

            for (auto & j: diffmodel->data_cols)
                exposure.at(i) += agent(j) * diffmodel->params.at(j);

            // Baseline probability of adoption
            double p = m->get_viruses()[i]->get_prob_infecting(m);
            exposure.at(i) += std::log(p) - std::log(1.0 - p);

            // Computing as log
            exposure.at(i) = 1.0/(1.0 + std::exp(-exposure.at(i)));

        }

        // Running the roulette to see is an innovation is adopted
        int which = roulette<int>(exposure, m);

        // No innovation was adopted
        if (which < 0)
            return;

        // Otherwise, it is adopted from any of the neighbors
        agent.set_virus(
            *innovations.at(which),
            m,
            ModelDiffNet::ADOPTER
        );

        return;

        };

    // Adding agents data
    model.set_agents_data(agents_data, data_ncols);
    
    // Adding statuses
    model.add_state("Non adopters", update_non_adopters);
    model.add_state("Adopters");

    // Adding parameters
    std::string parname = std::string("Prob. Adopting ") + innovation_name;
    model.add_param(prob_adopt, parname);

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> innovation(innovation_name, prevalence, true);
    innovation.set_state(1,1,1);
    
    innovation.set_prob_infecting(&model(parname));
    
    model.add_virus(innovation);

    model.set_name(
        std::string("Diffusion of Innovations - ") + innovation_name);

    return;
   
}

template<typename TSeq>
inline ModelDiffNet<TSeq>::ModelDiffNet(
    const std::string & innovation_name,
    epiworld_double prevalence,
    epiworld_double prob_adopt,
    bool normalize_exposure,
    double * agents_data,
    size_t data_ncols,
    std::vector< size_t > data_cols,
    std::vector< double > params
    )
{

    ModelDiffNet<TSeq>(
        *this,
        innovation_name,
        prevalence,
        prob_adopt,
        normalize_exposure,
        agents_data,
        data_ncols,
        data_cols,
        params
        );

    return;

}

#endif
