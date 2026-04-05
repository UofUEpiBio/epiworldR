#ifndef EPIWORLD_MODELS_SIRDCONNECTED_HPP 
#define EPIWORLD_MODELS_SIRDCONNECTED_HPP

#include "../model-bones.hpp"

/**
 * @brief Template for a Susceptible-Infected-Removed-Deceased (SIRD) model with connected population
 * 
 * ![Model Diagram](../assets/img/sirdconnected.png)
 * 
 * @ingroup connected_models
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRDCONN : public Model<TSeq>
{
public:
    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;
    static const int DECEASED    = 3;

    ModelSIRDCONN() {
        
        // tracked_agents_infected.reserve(1e4);
        // tracked_agents_infected_next.reserve(1e4);

    };

    ModelSIRDCONN(
        ModelSIRDCONN<TSeq> & model,
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate, 
        epiworld_double death_rate
    );

    ModelSIRDCONN(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );

    std::unique_ptr< Model<TSeq> > clone_ptr() override;


};

template<typename TSeq>
inline std::unique_ptr<Model<TSeq>> ModelSIRDCONN<TSeq>::clone_ptr()
{
    
    return std::make_unique<ModelSIRDCONN<TSeq>>(*this);

}

/**
 * @brief Template for a Susceptible-Infected-Removed-Deceased (SIRD) model
 * 
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 * @param death_rate Probability of death
 */
template<typename TSeq>
inline ModelSIRDCONN<TSeq>::ModelSIRDCONN(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    // epiworld_double prob_reinfection
    )
{

    UpdateFun<TSeq> update_susceptible = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void
        {

            // Sampling how many individuals
            m->set_rand_binom(
                m->size(),
                static_cast<double>(
                    m->par("Contact rate"))/
                    static_cast<double>(m->size())
            );

            int ndraw = m->rbinom();

            if (ndraw == 0)
                return;

            // Drawing from the set
            int nviruses_tmp = 0;
            auto & m_ref = *m;
            for (int i = 0; i < ndraw; ++i)
            {
                // Now selecting who is transmitting the disease
                int which = m->runif_int(0, static_cast<int>(m->size()) - 1);

                // Can't sample itself
                if (which == static_cast<int>(p->get_id()))
                    continue;

                // If the neighbor is infected, then proceed
                auto & neighbor = m->get_agents()[which];
                if (neighbor.get_state() == ModelSIRDCONN<TSeq>::INFECTED)
                {

                    auto & v = neighbor.get_virus();
                    
                    #ifdef EPI_DEBUG
                    if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
                        throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                    #endif
                        
                    /* And it is a function of susceptibility_reduction as well */ 
                    m->array_double_tmp[nviruses_tmp] =
                        (1.0 - p->get_susceptibility_reduction(v, m_ref)) *
                        v->get_prob_infecting(m) *
                        (1.0 - neighbor.get_transmission_reduction(v, m_ref))
                        ;
                
                    m->array_virus_tmp[nviruses_tmp++] = &(*v);

                }
            }

            // No virus to compute
            if (nviruses_tmp == 0u)
                return;

            // Running the roulette
            int which = roulette(nviruses_tmp, m);

            if (which < 0)
                return;

            p->set_virus(*m, *m->array_virus_tmp[which]);

            return; 

        };


    UpdateFun<TSeq> update_infected = [](
        Agent<TSeq> * p, Model<TSeq> * m
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSIRDCONN<TSeq>::INFECTED)
            {


                // Odd: Die, Even: Recover
                epiworld_fast_uint n_events = 0u;
                auto & v = p->get_virus();
                    
                // Die
                m->array_double_tmp[n_events++] = 
                v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, *m));
                
                // Recover
                m->array_double_tmp[n_events++] = 
                1.0 - (1.0 - v->get_prob_recovery(m)) *
                    (1.0 - p->get_recovery_enhancer(v, *m));
                
    #ifdef EPI_DEBUG
                if (n_events == 0u)
                {
                    printf_epiworld(
                    "[epi-debug] agent %i has 0 possible events!!\n",
                    static_cast<int>(p->get_id())
                    );
                    throw std::logic_error("Zero events in exposed.");
                }
    #else
                if (n_events == 0u)
                    return;
    #endif
                
                
                // Running the roulette
                int which = roulette(n_events, m);
                
                if (which < 0)
                    return;
                
                // Which roulette happen?
                if ((which % 2) == 0) // If odd
                {
                    
                    p->rm_virus(*m, ModelSIRDCONN<TSeq>::DECEASED);
                    
                } else {
                    
                    p->rm_virus(*m);
                    
                }

                return ;

            } else
                throw std::logic_error("This function can only be applied to infected individuals. (SIR)") ;

            return;

        };

    // state
    this->add_state("Susceptible", update_susceptible);
    this->add_state("Infected", update_infected);
    this->add_state("Recovered");
    this->add_state("Deceased");
      

    // Setting up parameters
    this->add_param(contact_rate, "Contact rate");
    this->add_param(transmission_rate, "Transmission rate");
    this->add_param(recovery_rate, "Recovery rate");
    this->add_param(death_rate, "Death rate");
    // this->add_param(prob_reinfection, "Prob. Reinfection");
    
    // Preparing the virus -------------------------------------------
    Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(1, 2, 3);
    virus.set_prob_infecting("Transmission rate");
    virus.set_prob_recovery("Recovery rate");
    virus.set_prob_death("Death rate");
    
    this->add_virus(virus);

    this->queuing_off(); // No queuing need

    this->agents_empty_graph(n);

    this->set_name("Susceptible-Infected-Removed-Deceased (SIRD) (connected)");

}


#endif
