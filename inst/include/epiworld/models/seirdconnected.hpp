#ifndef EPIWORLD_MODELS_SEIRDCONNECTED_HPP
#define EPIWORLD_MODELS_SEIRDCONNECTED_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRDCONN : public epiworld::Model<TSeq> 
{
public:

    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int REMOVED     = 3;
    static const int DECEASED    = 4;


    ModelSEIRDCONN() {};

    ModelSEIRDCONN(
        ModelSEIRDCONN<TSeq> & model,
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );
    
    ModelSEIRDCONN(
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        epiworld_double death_rate
    );

    void run(
        epiworld_fast_uint ndays,
        int seed = -1
    );

    void reset();

    Model<TSeq> * clone_ptr();

};

template<typename TSeq>
inline void ModelSEIRDCONN<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{
    
    Model<TSeq>::run(ndays, seed);

}

template<typename TSeq>
inline void ModelSEIRDCONN<TSeq>::reset()
{

    Model<TSeq>::reset();

    Model<TSeq>::set_rand_binom(
        Model<TSeq>::size(),
        static_cast<double>(
            Model<TSeq>::par("Contact rate"))/
            static_cast<double>(Model<TSeq>::size())
        );

    return;

}

template<typename TSeq>
inline Model<TSeq> * ModelSEIRDCONN<TSeq>::clone_ptr()
{
    
    ModelSEIRDCONN<TSeq> * ptr = new ModelSEIRDCONN<TSeq>(
        *dynamic_cast<const ModelSEIRDCONN<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 * @param death_rate Probability of death
 */
template<typename TSeq>
inline ModelSEIRDCONN<TSeq>::ModelSEIRDCONN(
    ModelSEIRDCONN<TSeq> & model,
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    // epiworld_double prob_reinfection
    )
{

    epiworld::UpdateFun<TSeq> update_susceptible = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void
        {

            // Sampling how many individuals
            int ndraw = m->rbinom();

            if (ndraw == 0)
                return;

            // Drawing from the set
            int nviruses_tmp = 0;
            for (int i = 0; i < ndraw; ++i)
            {
                // Now selecting who is transmitting the disease
                int which = static_cast<int>(
                    std::floor(m->size() * m->runif())
                );

                /* There is a bug in which runif() returns 1.0. It is rare, but
                 * we saw it here. See the Notes section in the C++ manual
                 * https://en.cppreference.com/mwiki/index.php?title=cpp/numeric/random/uniform_real_distribution&oldid=133329
                 * And the reported bug in GCC:
                 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
                 * 
                 */
                if (which == static_cast<int>(m->size()))
                    --which;

                // Can't sample itself
                if (which == static_cast<int>(p->get_id()))
                    continue;

                // If the neighbor is infected, then proceed
                auto & neighbor = m->get_agents()[which];
                if (neighbor.get_state() == ModelSEIRDCONN<TSeq>::INFECTED)
                {

                    for (const VirusPtr<TSeq> & v : neighbor.get_viruses()) 
                    { 

                        #ifdef EPI_DEBUG
                        if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
                            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                        #endif
                            
                        /* And it is a function of susceptibility_reduction as well */ 
                        m->array_double_tmp[nviruses_tmp] =
                            (1.0 - p->get_susceptibility_reduction(v, m)) * 
                            v->get_prob_infecting(m) * 
                            (1.0 - neighbor.get_transmission_reduction(v, m)) 
                            ; 
                    
                        m->array_virus_tmp[nviruses_tmp++] = &(*v);
                        
                    } 

                }
            }

            // No virus to compute
            if (nviruses_tmp == 0u)
                return;

            // Running the roulette
            int which = roulette(nviruses_tmp, m);

            if (which < 0)
                return;

            p->add_virus(
                *m->array_virus_tmp[which],
                m,
                ModelSEIRDCONN<TSeq>::EXPOSED
                );

            return; 

        };

    epiworld::UpdateFun<TSeq> update_infected = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSEIRDCONN<TSeq>::EXPOSED)
            {

                // Getting the virus
                auto & v = p->get_virus(0u);

                // Does the agent become infected?
                if (m->runif() < 1.0/(v->get_incubation(m)))
                {

                    p->change_state(m, ModelSEIRDCONN<TSeq>::INFECTED);
                    return;

                }


            } else if (state == ModelSEIRDCONN<TSeq>::INFECTED)
            {


              // Odd: Die, Even: Recover
              epiworld_fast_uint n_events = 0u;
              for (const auto & v : p->get_viruses())
              {
                
                // Die
                m->array_double_tmp[n_events++] = 
                  v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, m)); 
                
                // Recover
                m->array_double_tmp[n_events++] = 
                  1.0 - (1.0 - v->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(v, m)); 
                
              }
              
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
                
                size_t which_v = std::ceil(which / 2);
                p->rm_agent_by_virus(which_v, m);
                
              } else {
                
                size_t which_v = std::floor(which / 2);
                p->rm_virus(which_v, m);
                
              }

                return ;

            } else
                throw std::logic_error("This function can only be applied to exposed or infected individuals. (SEIRD)") ;

            return;

        };

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Prob. Transmission");
    model.add_param(recovery_rate, "Prob. Recovery");
    model.add_param(avg_incubation_days, "Avg. Incubation days");
    model.add_param(death_rate, "Death rate");
    
    // state
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Exposed", update_infected);
    model.add_state("Infected", update_infected);
    model.add_state("Removed");
    model.add_state("Deceased");


    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(
        ModelSEIRDCONN<TSeq>::EXPOSED,
        ModelSEIRDCONN<TSeq>::REMOVED,
        ModelSEIRDCONN<TSeq>::DECEASED
        );

    virus.set_prob_infecting(&model("Prob. Transmission"));
    virus.set_prob_recovery(&model("Prob. Recovery"));
    virus.set_incubation(&model("Avg. Incubation days"));
    virus.set_prob_death(&model("Death rate"));
    model.add_virus(virus, prevalence);

    model.queuing_off(); // No queuing need

    // Adding the empty population
    model.agents_empty_graph(n);

    model.set_name("Susceptible-Exposed-Infected-Removed-Deceased (SEIRD) (connected)");

    return;

}

template<typename TSeq>
inline ModelSEIRDCONN<TSeq>::ModelSEIRDCONN(
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    epiworld_double death_rate
    )
{

    ModelSEIRDCONN(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        transmission_rate,
        avg_incubation_days,
        recovery_rate,
        death_rate
    );

    return;

}

#endif
