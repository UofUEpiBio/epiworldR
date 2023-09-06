#ifndef EPIWORLD_MODELS_SEIRCONNECTEDLOGIT_HPP
#define EPIWORLD_MODELS_SEIRCONNECTEDLOGIT_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSEIRCONNLogit : public epiworld::Model<TSeq> 
{
private:
    static const int SUSCEPTIBLE = 0;
    static const int EXPOSED     = 1;
    static const int INFECTED    = 2;
    static const int RECOVERE    = 3;

public:

    ModelSEIRCONNLogit() {
        tracked_agents_infected.reserve(1e4);
        tracked_agents_infected_next.reserve(1e4);
    };

    ModelSEIRCONNLogit(
        ModelSEIRCONNLogit<TSeq> & model,
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate,
        double * covars,
        std::vector< double > logit_params
    );
    
    ModelSEIRCONNLogit(
        std::string vname,
        epiworld_fast_uint n,
        epiworld_double prevalence,
        epiworld_double contact_rate,
        epiworld_double transmission_rate,
        epiworld_double avg_incubation_days,
        epiworld_double recovery_rate
        double * covars,
        std::vector< double > logit_params
    );

    // Tracking who is infected and who is not
    std::vector< epiworld::Agent<>* > tracked_agents_infected = {};
    std::vector< epiworld::Agent<>* > tracked_agents_infected_next = {};

    bool tracked_started = false;
    int tracked_ninfected = 0;
    int tracked_ninfected_next = 0;

};

/**
 * @brief Template for a Susceptible-Exposed-Infected-Removed (SEIR) model
 * 
 * @param model A Model<TSeq> object where to set up the SIR.
 * @param vname std::string Name of the virus
 * @param prevalence Initial prevalence (proportion)
 * @param contact_rate Reproductive number (beta)
 * @param transmission_rate Probability of transmission
 * @param recovery_rate Probability of recovery
 */
template<typename TSeq>
inline ModelSEIRCONNLogit<TSeq>::ModelSEIRCONNLogit(
    ModelSEIRCONNLogit<TSeq> & model,
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    double * covars,
    std::vector< double > logit_params
    // epiworld_double prob_reinfection
    )
{

    auto * _tracked_agents_infected = &model.tracked_agents_infected;
    auto * _tracked_agents_infected_next = &model.tracked_agents_infected_next;
    auto * _tracked_started = &model.tracked_started;
    auto * _tracked_ninfected = &model.tracked_ninfected;
    auto * _tracked_ninfected_next = &model.tracked_ninfected_next;

    std::function<void(epiworld::Model<TSeq> *)> tracked_agents_check_init = 
    [
        _tracked_started,
        _tracked_agents_infected,
        _tracked_ninfected
    ](epiworld::Model<TSeq> * m) 
        {

            if (*_tracked_started)
                return;

            /* Checking first if it hasn't  */ 
            if (!*_tracked_started) 
            { 
                
                /* Listing who is infected */ 
                for (auto & p : m->get_agents())
                {
                    if (p.get_state() == ModelSEIRCONNLogit<TSeq>::INFECTED)
                    {
                    
                        _tracked_agents_infected->push_back(&p);
                        *_tracked_ninfected += 1;
                    
                    }
                }

                for (auto & p: *_tracked_agents_infected)
                {
                    if (p->get_n_viruses() == 0)
                        throw std::logic_error("Cannot be infected and have no viruses.");
                }
                
                *_tracked_started = true;
                
            }

        };

    epiworld::UpdateFun<TSeq> update_susceptible = 
    [
        tracked_agents_check_init,
        _tracked_ninfected,
        _tracked_agents_infected
    ](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
        {

            tracked_agents_check_init(m);

            // No infected individual?
            if (*_tracked_ninfected == 0)
                return;

            // Computing probability of contagion
            // P(infected) = 1 - (1 - beta/Pop * ptransmit) ^ ninfected
            epiworld_double prob_infect = 1.0 - std::pow(
                1.0 - (m->par("Contact rate")) * (m->par("Transmission rate")) / m->size(),
                *_tracked_ninfected
                );

            if (m->runif() < prob_infect)
            {

                // Now selecting who is transmitting the disease
                epiworld_fast_uint which = static_cast<epiworld_fast_uint>(
                    std::floor(*_tracked_ninfected * m->runif())
                );

                // Infecting the individual
                #ifdef EPI_DEBUG
                if (_tracked_agents_infected->operator[](which)->get_n_viruses() == 0)
                {

                    printf_epiworld("[epi-debug] date: %i\n", m->today());
                    printf_epiworld("[epi-debug] sim#: %i\n", m->get_n_replicates());

                    throw std::logic_error(
                        "[epi-debug] The agent " + std::to_string(which) + " has no "+
                        "virus to share. The agent's state is: " +
                        std::to_string(_tracked_agents_infected->operator[](which)->get_state())
                    );
                }
                #endif
                p->add_virus(
                    _tracked_agents_infected->operator[](which)->get_virus(0u),
                    ModelSEIRCONNLogit<TSeq>::EXPOSED
                    ); 

                return;

            }

            return;

        };

    epiworld::UpdateFun<TSeq> update_infected = 
    [
        tracked_agents_check_init,
        _tracked_agents_infected_next,
        _tracked_ninfected_next

    ](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
        {

            tracked_agents_check_init(m);
            auto state = p->get_state();

            if (state == ModelSEIRCONNLogit<TSeq>::EXPOSED)
            {

                // Does the agent become infected?
                if (m->runif() < 1.0/(m->par("Avg. Incubation days")))
                {
                    // Adding the individual to the queue
                    _tracked_agents_infected_next->push_back(p);
                    *_tracked_ninfected_next += 1;

                    p->change_state(m, ModelSEIRCONNLogit<TSeq>::INFECTED);

                    return;

                }


            } else if (state == ModelSEIRCONNLogit<TSeq>::INFECTED)
            {

                if (m->runif() < (m->par("Recovery rate")))
                {

                    *_tracked_ninfected_next -= 1;
                    p->rm_virus(0, m);
                    return;

                }

                _tracked_agents_infected_next->push_back(p);

            } 

            return;

        };

    epiworld::GlobalFun<TSeq> global_accounting = 
    [
        _tracked_started,
        _tracked_agents_infected,
        _tracked_agents_infected_next,
        _tracked_ninfected,
        _tracked_ninfected_next
    ](epiworld::Model<TSeq>* m) -> void
        {

            // On the last day, also reset tracked agents and
            // set the initialized value to false
            if (static_cast<epiworld_fast_uint>(m->today()) == (m->get_ndays() - 1))
            {

                *_tracked_started = false;
                _tracked_agents_infected->clear();
                _tracked_agents_infected_next->clear();
                *_tracked_ninfected = 0;
                *_tracked_ninfected_next = 0;    

                return;
            }

            std::swap(*_tracked_agents_infected, *_tracked_agents_infected_next);
            _tracked_agents_infected_next->clear();

            *_tracked_ninfected += *_tracked_ninfected_next;
            *_tracked_ninfected_next = 0;

        };

    // Setting up parameters
    model.add_param(contact_rate, "Contact rate");
    model.add_param(transmission_rate, "Transmission rate");
    model.add_param(recovery_rate, "Recovery rate");
    model.add_param(avg_incubation_days, "Avg. Incubation days");
    
    // state
    model.add_state("Susceptible", update_susceptible);
    model.add_state("Exposed", update_infected);
    model.add_state("Infected", update_infected);
    model.add_state("Recovered");

    // Adding agent's parameters
    model.set_agents_data(covars, logit_params.size());
    for (size_t i = 0u; i < logit_params.size(); ++i)
        model.add_param(
            std::string("x" + std::to_string(i)),
            logit_params[i]
        );

    // Preparing the virus -------------------------------------------
    epiworld::Virus<TSeq> virus(vname);
    virus.set_state(1,3,3);
    model.add_virus(virus, prevalence);

    // Adding updating function
    model.add_global_action(global_accounting, "Accounting", -1);

    model.queuing_off(); // No queuing need

    // Adding the empty population
    model.agents_empty_graph(n);

    model.set_name("Susceptible-Exposed-Infected-Removed (SEIR) (connected)");

    return;

}

template<typename TSeq>
inline ModelSEIRCONNLogit<TSeq>::ModelSEIRCONNLogit(
    std::string vname,
    epiworld_fast_uint n,
    epiworld_double prevalence,
    epiworld_double contact_rate,
    epiworld_double transmission_rate,
    epiworld_double avg_incubation_days,
    epiworld_double recovery_rate,
    double * covars,
    std::vector< double > logit_params
    )
{

    ModelSEIRCONNLogit(
        *this,
        vname,
        n,
        prevalence,
        contact_rate,
        transmission_rate,
        avg_incubation_days,
        covars,
        logit_params
    );

    return;

}

#endif
