#ifndef EPIWORLD_MODELS_SURVEILLANCE_HPP
#define EPIWORLD_MODELS_SURVEILLANCE_HPP

template<typename TSeq = EPI_DEFAULT_TSEQ>
class ModelSURV : public epiworld::Model<TSeq> {

private:
    // state
    static const int SUSCEPTIBLE           = 0;
    static const int LATENT                = 1;
    static const int SYMPTOMATIC           = 2;
    static const int SYMPTOMATIC_ISOLATED  = 3; // sampled and discovered
    static const int ASYMPTOMATIC          = 4;
    static const int ASYMPTOMATIC_ISOLATED = 5;
    static const int RECOVERED             = 6;
    static const int REMOVED               = 7;

public:

    /**
     * @brief Vector of days spent in latent and infectious states
     * A row-major matrix
     */
    std::vector< epiworld_double > days_latent_and_infectious;

    /**
     * @name Construct a new ModelSURV object
     * 
     * The ModelSURV class simulates a survaillence model where agents can be
     * isolated, even if asyptomatic.
     * 
     * @param vname String. Name of the virus
     * @param prevalence Integer. Number of initial cases of the virus.
     * @param efficacy_vax Double. Efficacy of the vaccine (1 - P(acquire the disease)).
     * @param latent_period Double. Shape parameter of a `Gamma(latent_period, 1)`
     *   distribution. This coincides with the expected number of latent days.
     * @param infect_period Double. Shape parameter of a `Gamma(infected_period, 1)`
     *   distribution. This coincides with the expected number of infectious days.
     * @param prob_symptoms Double. Probability of generating symptoms.
     * @param prop_vaccinated Double. Probability of vaccination. Coincides with
     *   the initial prevalence of vaccinated individuals.
     * @param prop_vax_redux_transm Double. Factor by which the vaccine reduces
     *   transmissibility.
     * @param prop_vax_redux_infect Double. Factor by which the vaccine reduces
     *   the chances of becoming infected.
     * @param surveillance_prob Double. Probability of testing an agent.
     * @param prob_transmission Double. Raw transmission probability.
     * @param prob_death Double. Raw probability of death for symptomatic individuals.
     * @param prob_noreinfect Double. Probability of no re-infection.
     * 
     * @details
     * This model features the following states:
     * 
     * - Susceptible
     * - Latent
     * - Symptomatic
     * - Symptomatic isolated
     * - Asymptomatic
     * - Asymptomatic isolated
     * - Recovered
     * - Removed    
     * 
     * @returns An object of class `epiworld_surv`
     * 
     */
    ///@{
    ModelSURV() {};

    ModelSURV(
        ModelSURV<TSeq> & model,
        const std::string & vname,
        epiworld_fast_uint prevalence               = 50,
        epiworld_double efficacy_vax          = 0.9,
        epiworld_double latent_period         = 3u,
        epiworld_double infect_period         = 6u,
        epiworld_double prob_symptoms         = 0.6,
        epiworld_double prop_vaccinated       = 0.25,
        epiworld_double prop_vax_redux_transm = 0.5,
        epiworld_double prop_vax_redux_infect = 0.5,
        epiworld_double surveillance_prob     = 0.001,
        epiworld_double prob_transmission     = 1.0,
        epiworld_double prob_death            = 0.001,
        epiworld_double prob_noreinfect       = 0.9
    );

    ModelSURV(
        const std::string & vname,
        epiworld_fast_uint prevalence         = 50,
        epiworld_double efficacy_vax          = 0.9,
        epiworld_double latent_period         = 3u,
        epiworld_double infect_period         = 6u,
        epiworld_double prob_symptoms         = 0.6,
        epiworld_double prop_vaccinated       = 0.25,
        epiworld_double prop_vax_redux_transm = 0.5,
        epiworld_double prop_vax_redux_infect = 0.5,
        epiworld_double surveillance_prob     = 0.001,
        epiworld_double prob_transmission     = 1.0,
        epiworld_double prob_death            = 0.001,
        epiworld_double prob_noreinfect       = 0.9
    );
    ///@}

    void reset();

};

template<typename TSeq>
inline void ModelSURV<TSeq>::reset()
{
    epiworld::Model<TSeq>::reset();

    days_latent_and_infectious.clear();
    days_latent_and_infectious.resize(
        2u * epiworld::Model<TSeq>::size(),
        -1.0
    );
}

template<typename TSeq>
inline ModelSURV<TSeq>::ModelSURV(
    ModelSURV<TSeq> & model,
    const std::string & vname,
    epiworld_fast_uint prevalence,
    epiworld_double efficacy_vax,
    epiworld_double latent_period,
    epiworld_double infect_period,
    epiworld_double prob_symptoms,
    epiworld_double prop_vaccinated,
    epiworld_double prop_vax_redux_transm,
    epiworld_double prop_vax_redux_infect,
    epiworld_double surveillance_prob,
    epiworld_double prob_transmission,
    epiworld_double prob_death,
    epiworld_double prob_noreinfect
    )
{

    EPI_NEW_UPDATEFUN_LAMBDA(surveillance_update_susceptible, TSeq) {

        // This computes the prob of getting any neighbor variant
        epiworld_fast_uint nviruses_tmp = 0u;
        for (auto & neighbor: p->get_neighbors()) 
        {
                    
            auto & v = neighbor->get_virus();

            if (v == nullptr)
                continue;
                
            /* And it is a function of susceptibility_reduction as well */ 
            epiworld_double tmp_transmission = 
                (1.0 - p->get_susceptibility_reduction(v, m)) * 
                v->get_prob_infecting(m) * 
                (1.0 - neighbor->get_transmission_reduction(v, m)) 
                ; 
        
            m->array_double_tmp[nviruses_tmp]  = tmp_transmission;
            m->array_virus_tmp[nviruses_tmp++] = &(*v);
        }

        // No virus to compute on
        if (nviruses_tmp == 0)
            return;

        // Running the roulette
        int which = roulette(nviruses_tmp, m);

        if (which < 0)
            return;

        p->set_virus(*m->array_virus_tmp[which], m); 
        return;

    };


    epiworld::UpdateFun<TSeq> surveillance_update_exposed = 
    [](epiworld::Agent<TSeq> * p, epiworld::Model<TSeq> * m) -> void
    {

        // Dynamically getting the ModelSURV
        ModelSURV<TSeq> * model_surv = dynamic_cast<ModelSURV<TSeq> *>(m);

        epiworld::VirusPtr<TSeq> & v = p->get_virus(); 
        epiworld_double p_die = v->get_prob_death(m) * (1.0 - p->get_death_reduction(v, m)); 
        
        epiworld_fast_uint days_since_exposed = m->today() - v->get_date();
        epiworld_fast_uint state = p->get_state();

        // Figuring out latent period
        auto & dat = model_surv->days_latent_and_infectious;
        if (dat[p->get_id()] < 0)
        {
            epiworld_double latent_days = m->rgamma(
                m->par("Latent period"), 1.0
            );

            dat[p->get_id() * 2u] = latent_days;

            dat[p->get_id() * 2u + 1u] = 
                m->rgamma(m->par("Infect period"), 1.0) +
                latent_days;
        }
        
        // If still latent, nothing happens
        if (days_since_exposed <= dat[p->get_id() * 2u])
            return;

        // If past days infected + latent, then bye.
        if (days_since_exposed >= dat[p->get_id() * 2u + 1u])
        {
            p->rm_virus(m);
            return;
        }

        // If it is infected, then it can be asymptomatic or symptomatic
        if (state == ModelSURV<TSeq>::LATENT)
        {

            // Will be symptomatic?
            if (EPI_RUNIF() < m->par("Prob of symptoms"))
                p->change_state(m, ModelSURV<TSeq>::SYMPTOMATIC);
            else
                p->change_state(m, ModelSURV<TSeq>::ASYMPTOMATIC);
            
            return;

        }
        
        // Otherwise, it can be removed
        if (EPI_RUNIF() < p_die)
        {
            p->change_state(m, ModelSURV<TSeq>::REMOVED, -1);
            return;
        }
        
        return;

    };

    std::vector< unsigned int > exposed_state = {
        SYMPTOMATIC,
        SYMPTOMATIC_ISOLATED,
        ASYMPTOMATIC,
        ASYMPTOMATIC_ISOLATED,
        LATENT
    };

    epiworld::GlobalFun<TSeq> surveillance_program = 
    [exposed_state](
        epiworld::Model<TSeq>* m
        ) -> void
    {

        // How many will we find
        std::binomial_distribution<> bdist(m->size(), m->par("Surveilance prob."));
        int nsampled = bdist(*m->get_rand_endgine());

        int to_go = nsampled + 1;

        epiworld_double ndetected        = 0.0;
        epiworld_double ndetected_asympt = 0.0;
        
        auto & pop = m->get_agents();
        std::vector< bool > sampled(m->size(), false);
        
        while (to_go-- > 0)
        {

            // Who is the lucky one
            epiworld_fast_uint i = static_cast<epiworld_fast_uint>(std::floor(EPI_RUNIF() * m->size()));

            if (sampled[i])
                continue;

            sampled[i] = true;
            epiworld::Agent<TSeq> * p = &pop[i];
            
            // If still exposed for the next term
            if (epiworld::IN(p->get_state(), exposed_state ))
            {

                ndetected += 1.0;
                if (p->get_state() == ModelSURV<TSeq>::ASYMPTOMATIC)
                {
                    ndetected_asympt += 1.0;
                    p->change_state(m, ModelSURV<TSeq>::ASYMPTOMATIC_ISOLATED);
                }
                else 
                {
                    p->change_state(m, ModelSURV<TSeq>::SYMPTOMATIC_ISOLATED);
                }

            }

        }

        // Writing the user data
        std::vector< int > totals;
        m->get_db().get_today_total(nullptr,&totals);
        m->add_user_data(
            {
                static_cast<epiworld_double>(nsampled),
                ndetected,
                ndetected_asympt,
                static_cast<epiworld_double>(totals[ModelSURV<TSeq>::ASYMPTOMATIC])
            }
            );


    };

    model.add_state("Susceptible", surveillance_update_susceptible);
    model.add_state("Latent", surveillance_update_exposed);
    model.add_state("Symptomatic", surveillance_update_exposed);
    model.add_state("Symptomatic isolated", surveillance_update_exposed);
    model.add_state("Asymptomatic", surveillance_update_exposed);
    model.add_state("Asymptomatic isolated", surveillance_update_exposed);
    model.add_state("Recovered");
    model.add_state("Removed");

    // General model parameters
    model.add_param(latent_period, "Latent period");
    model.add_param(infect_period, "Infect period");
    model.add_param(prob_symptoms, "Prob of symptoms");
    model.add_param(surveillance_prob, "Surveilance prob.");
    model.add_param(efficacy_vax, "Vax efficacy");
    model.add_param(prop_vax_redux_transm, "Vax redux transmission");
    model.add_param(prob_transmission, "Prob of transmission");
    model.add_param(prob_death, "Prob. death");
    model.add_param(prob_noreinfect, "Prob. no reinfect");

    // Virus ------------------------------------------------------------------
    epiworld::Virus<TSeq> covid("Covid19", prevalence, false);
    covid.set_state(LATENT, RECOVERED, REMOVED);
    covid.set_post_immunity(&model("Prob. no reinfect"));
    covid.set_prob_death(&model("Prob. death"));

    epiworld::VirusFun<TSeq> ptransmitfun = [](
        epiworld::Agent<TSeq> * p,
        epiworld::Virus<TSeq> &,
        epiworld::Model<TSeq> * m
        ) -> epiworld_double
    {
        // No chance of infecting
        epiworld_fast_uint  s = p->get_state();
        if (s == ModelSURV<TSeq>::LATENT)
            return static_cast<epiworld_double>(0.0);
        else if (s == ModelSURV<TSeq>::SYMPTOMATIC_ISOLATED)
            return static_cast<epiworld_double>(0.0);
        else if (s == ModelSURV<TSeq>::ASYMPTOMATIC_ISOLATED)
            return static_cast<epiworld_double>(0.0);

        // Otherwise
        return m->par("Prob of transmission");
    };

    covid.set_prob_infecting_fun(ptransmitfun);
    
    model.add_virus(covid);

    model.set_user_data({"nsampled", "ndetected", "ndetected_asympt", "nasymptomatic"});
    model.add_globalevent(surveillance_program, "Surveilance program", -1);
   
    // Vaccine tool -----------------------------------------------------------
    epiworld::Tool<TSeq> vax("Vaccine", prop_vaccinated, true);
    vax.set_susceptibility_reduction(&model("Vax efficacy"));
    vax.set_transmission_reduction(&model("Vax redux transmission"));
    
    model.add_tool(vax);

    model.set_name("Surveillance");

    return;

}

template<typename TSeq>
inline ModelSURV<TSeq>::ModelSURV(
    const std::string & vname,
    epiworld_fast_uint prevalence,
    epiworld_double efficacy_vax,
    epiworld_double latent_period,
    epiworld_double infect_period,
    epiworld_double prob_symptoms,
    epiworld_double prop_vaccinated,
    epiworld_double prop_vax_redux_transm,
    epiworld_double prop_vax_redux_infect,
    epiworld_double surveillance_prob,
    epiworld_double prob_transmission,
    epiworld_double prob_death,
    epiworld_double prob_noreinfect
    )
{

    ModelSURV(
        *this,
        vname,
        prevalence,
        efficacy_vax,
        latent_period,
        infect_period,
        prob_symptoms,
        prop_vaccinated,
        prop_vax_redux_transm,
        prop_vax_redux_infect,
        surveillance_prob,
        prob_transmission,
        prob_death,
        prob_noreinfect
    );

    return;

}

#endif
