#ifndef EPIWORLD_MODELS_SIRCONNECTED_HPP 
#define EPIWORLD_MODELS_SIRCONNECTED_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSIRCONN : public epiworld::Model<TSeq>
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int INFECTED    = 1;
    static const int RECOVERED   = 2;

public:

    ModelSIRCONN() {};

    ModelSIRCONN(
        ModelSIRCONN<TSeq> & model,
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate
    );

    ModelSIRCONN(
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double recovery_rate
    );

    ModelSIRCONN<TSeq> & run(
        epiworld_fast_uint ndays,
        int seed = -1
    );
    
    void reset();

    Model<TSeq> * clone_ptr();

    /**
     * @brief Set the initial states of the model
     * @param proportions_ Double vector with a single element:
     * - The proportion of non-infected individuals who have recovered.
    */
    ModelSIRCONN<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    );


};

template<typename TSeq>
inline ModelSIRCONN<TSeq> & ModelSIRCONN<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{

    Model<TSeq>::run(ndays, seed);
    return *this;

}

template<typename TSeq>
inline void ModelSIRCONN<TSeq>::reset()
{

    Model<TSeq>::reset();
    return;

}

template<typename TSeq>
inline Model<TSeq> * ModelSIRCONN<TSeq>::clone_ptr()
{
    
    ModelSIRCONN<TSeq> * ptr = new ModelSIRCONN<TSeq>(
        *dynamic_cast<const ModelSIRCONN<TSeq>*>(this)
        );

    return dynamic_cast< Model<TSeq> *>(ptr);

}

/**
 * @brief Template for a Susceptible-Infected-Removed (SIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Average number of contacts (interactions) per step.
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 */
template<typename TSeq>
inline ModelSIRCONN<TSeq>::ModelSIRCONN(
    ModelSIRCONN<TSeq> & model,
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate
    // epiworld_double prob_reinfection
    )
{

    epiworld::UpdateFun<TSeq> update_susceptible = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
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
                if (neighbor.get_state() == ModelSIRCONN<TSeq>::INFECTED)
                {

                    if (neighbor.get_virus() == nullptr)
                        continue;

                    auto & v = neighbor.get_virus();

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
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSIRCONN<TSeq>::INFECTED)
            {


                // Odd: Die, Even: Recover
                epiworld_fast_uint n_events = 0u;
                // Recover
                m->array_double_tmp[n_events++] = 
                    1.0 - (1.0 - p->get_virus()->get_prob_recovery(m)) *
                        (1.0 - p->get_recovery_enhancer(p->get_virus(), m)); 

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
                p->rm_virus(m);

                return ;

            } else
                throw std::logic_error(
                    "This function can only be applied to infected individuals. (SIR)"
                    ) ;

            return;

        };

    // state
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(recovery_rate, "Recovery rate");
    // model.add_param(prob_reinfection, "Prob. Reinfection");
    
    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1, 2, 2);
    virus.set_prob_infecting(&model("Transmission rate"));
    virus.set_prob_recovery(&model("Recovery rate"));

    model.add_virus(virus, prevalence);

    model.queuing_off(); // No queuing need

    model.agents_empty_graph(n);

    model.set_name("Susceptible-Infected-Removed (SIR) (connected)");

    return;

}

template<typename TSeq>
inline ModelSIRCONN<TSeq>::ModelSIRCONN(
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double recovery_rate
    )
{

    ModelSIRCONN(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        transmission_rate,
        recovery_rate
    );

    return;

}

template<typename TSeq>
inline ModelSIRCONN<TSeq> & ModelSIRCONN<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /**/ 
) {

    Model<TSeq>::initial_states_fun = 
        create_init_function_sir<TSeq>(proportions_)
        ;

    return *this;

}


#endif
