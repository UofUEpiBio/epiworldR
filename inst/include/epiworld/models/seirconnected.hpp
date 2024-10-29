#ifndef EPIWORLD_MODELS_SEIRCONNECTED_HPP
#define EPIWORLD_MODELS_SEIRCONNECTED_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRCONN : public epiworld::Model<TSeq> 
{
private:
    std::vector< epiworld::Agent<TSeq> * > infected;
    void update_infected();

public:

    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int RECOVERED   = 3;


    ModelSEIRCONN() {};

    ModelSEIRCONN(
        ModelSEIRCONN<TSeq> & model,
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate
    );
    
    ModelSEIRCONN(
        const std::string & vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate
    );

    ModelSEIRCONN<TSeq> & run(
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
    ModelSEIRCONN<TSeq> & initial_states(
        std::vector< double > proportions_,
        std::vector< int > queue_ = {}
    );

    size_t get_n_infected() const { return infected.size(); }

    /***
     * @brief Compute expected generation time
     * @param max_days Maximum number of days.
     * @param max_contacts Maximum number of contacts.
     */
    std::vector< double > generation_time_expected(
        int max_days = 200,
        int max_contacts = 200
    ) const;

};

template<typename TSeq>
inline void ModelSEIRCONN<TSeq>::update_infected()
{

    infected.clear();
    infected.reserve(this->size());

    for (auto & p : this->get_agents())
    {
        if (p.get_state() == ModelSEIRCONN<TSeq>::INFECTED)
        {
            infected.push_back(&p);
        }
    }

    Model<TSeq>::set_rand_binom(
        this->get_n_infected(),
        static_cast<double>(Model<TSeq>::par("Contact rate"))/
            static_cast<double>(Model<TSeq>::size())
    );

    return;

}

template<typename TSeq>
inline ModelSEIRCONN<TSeq> & ModelSEIRCONN<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
)
{
    
    Model<TSeq>::run(ndays, seed);

    return *this;

}

template<typename TSeq>
inline void ModelSEIRCONN<TSeq>::reset()
{

    Model<TSeq>::reset();
    this->update_infected();

    return;

}

template<typename TSeq>
inline Model<TSeq> * ModelSEIRCONN<TSeq>::clone_ptr()
{
    
    ModelSEIRCONN<TSeq> * ptr = new ModelSEIRCONN<TSeq>(
        *dynamic_cast<const ModelSEIRCONN<TSeq>*>(this)
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
 */
template<typename TSeq>
inline ModelSEIRCONN<TSeq>::ModelSEIRCONN(
    ModelSEIRCONN<TSeq> & model,
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate
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

            ModelSEIRCONN<TSeq> * model = dynamic_cast<ModelSEIRCONN<TSeq> *>(m);
            size_t ninfected = model->get_n_infected();

            // Drawing from the set
            int nviruses_tmp = 0;
            for (int i = 0; i < ndraw; ++i)
            {
                // Now selecting who is transmitting the disease
                int which = static_cast<int>(
                    std::floor(ninfected * m->runif())
                );

                /* There is a bug in which runif() returns 1.0. It is rare, but
                 * we saw it here. See the Notes section in the C++ manual
                 * https://en.cppreference.com/mwiki/index.php?title=cpp/numeric/random/uniform_real_distribution&oldid=133329
                 * And the reported bug in GCC:
                 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
                 * 
                 */
                if (which == static_cast<int>(ninfected))
                    --which;

                epiworld::Agent<TSeq> & neighbor = *model->infected[which];

                // Can't sample itself
                if (neighbor.get_id() == p->get_id())
                    continue;

                // The neighbor is infected by construction
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

            // No virus to compute
            if (nviruses_tmp == 0u)
                return;

            // Running the roulette
            int which = roulette(nviruses_tmp, m);

            if (which < 0)
                return;

            p->set_virus(
                *m->array_virus_tmp[which],
                m,
                ModelSEIRCONN<TSeq>::EXPOSED
                );

            return; 

        };

    epiworld::UpdateFun<TSeq> update_infected = [](
        epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m
        ) -> void {

            auto state = p->get_state();

            if (state == ModelSEIRCONN<TSeq>::EXPOSED)
            {

                // Getting the virus
                auto & v = p->get_virus();

                // Does the agent become infected?
                if (m->runif() < 1.0/(v->get_incubation(m)))
                {

                    p->change_state(m, ModelSEIRCONN<TSeq>::INFECTED);
                    return;

                }


            } else if (state == ModelSEIRCONN<TSeq>::INFECTED)
            {


                // Odd: Die, Even: Recover
                epiworld_fast_uint n_events = 0u;
                const auto & v = p->get_virus();

                // Recover
                m->array_double_tmp[n_events++] = 
                    1.0 - (1.0 - v->get_prob_recovery(m)) * (1.0 - p->get_recovery_enhancer(v, m)); 

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
                throw std::logic_error("This function can only be applied to exposed or infected individuals. (SEIR)") ;

            return;

        };

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Prob. Transmission");
    model.add_param(recovery_rate, "Prob. Recovery");
    model.add_param(avg_incubation_days, "Avg. Incubation days");
    
    // state
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Exposed", update_infected);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Adding update function
    epiworld::GlobalFun<TSeq> update = [](
        epiworld::Model<TSeq> * m
        ) -> void
        {

            ModelSEIRCONN<TSeq> * model = dynamic_cast<ModelSEIRCONN<TSeq> *>(m);

            model->update_infected();

            return;

        };

    model.add_globalevent(update, "Update infected individuals");


    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname, prevalence, true);
    virus.set_state(
        ModelSEIRCONN<TSeq>::EXPOSED,
        ModelSEIRCONN<TSeq>::RECOVERED,
        ModelSEIRCONN<TSeq>::RECOVERED
        );

    virus.set_prob_infecting(&model("Prob. Transmission"));
    virus.set_prob_recovery(&model("Prob. Recovery"));
    virus.set_incubation(&model("Avg. Incubation days"));

    model.add_virus(virus);

    model.queuing_off(); // No queuing need

    // Adding the empty population
    model.agents_empty_graph(n);

    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR) (connected)");

    return;

}

template<typename TSeq>
inline ModelSEIRCONN<TSeq>::ModelSEIRCONN(
    const std::string & vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate
    )
{

    ModelSEIRCONN(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        transmission_rate,
        avg_incubation_days,
        recovery_rate
    );

    return;

}

template<typename TSeq>
inline ModelSEIRCONN<TSeq> & ModelSEIRCONN<TSeq>::initial_states(
    std::vector< double > proportions_,
    std::vector< int > /* queue_ */
)
{

    Model<TSeq>::initial_states_fun =
        create_init_function_seir<TSeq>(proportions_)
        ;

    return *this;

}

template<typename TSeq>
inline std::vector< double > ModelSEIRCONN<TSeq>::generation_time_expected(
    int max_days,
    int max_contacts
) const
{

    // Retrieving total counts
    std::vector< int > h_date;
    std::vector< std::string > h_state;
    std::vector< int > h_counts;
    const auto this_const = dynamic_cast<const ModelSEIRCONN<TSeq> *>(this);
    this_const->get_db().get_hist_total(
        &h_date,
        &h_state,
        &h_counts
    );

    // Retrieving information on susceptibles
    std::vector< double > S(this_const->get_ndays(), 0.0);
    for (size_t i = 0; i < h_date.size(); ++i)
    {
        if (h_state[i] == "Susceptible")
            S[h_date[i]] += h_counts[i];
    }

    // Computing the expected number of days in exposed
    double days_exposed = this_const->par("Avg. Incubation days");

    // The generation time in the SEIR model starts from 2, as agents 
    // spend at least one day in the exposed state, and 1 day in the 
    // infectious state before starting transmitting.
    std::vector< double > gen_times(
        this_const->get_ndays(), 1.0 + days_exposed
        );
        
    double p_c = this_const->par("Contact rate")/this_const->size();
    double p_i = this_const->par("Prob. Transmission");
    double p_r = this_const->par("Prob. Recovery");

    for (size_t i = 0u; i < this_const->get_ndays(); ++i)
    {
        gen_times[i] += gen_int_mean(
            S[i],
            p_c,
            p_i,
            p_r,
            max_days,
            max_contacts
        );

    }

    return gen_times;

}

#endif
