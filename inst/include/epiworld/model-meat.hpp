#ifndef EPIWORLD_MODEL_MEAT_HPP
#define EPIWORLD_MODEL_MEAT_HPP

#include <vector>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <map>
#include <unordered_map>
#include "config.hpp"
#include "userdata-bones.hpp"
#include "adjlist-bones.hpp"
#include "model-bones.hpp"
#include "entities-bones.hpp"
#include "virus-bones.hpp"
#include "agent-bones.hpp"

/**
 * @brief Function factory for saving model runs
 * 
 * @details This function is the default behavior of the `run_multiple`
 * member of `Model<TSeq>`. By default only the total history (
 * case counts by unit of time.)
 * 
 * @tparam TSeq 
 * @param fmt 
 * @param total_hist 
 * @param virus_info 
 * @param virus_hist 
 * @param tool_info 
 * @param tool_hist 
 * @param transmission 
 * @param transition 
 * @return std::function<void(size_t,Model<TSeq>*)> 
 */
template<typename TSeq>
inline std::function<void(size_t,Model<TSeq>*)> make_save_run(
    std::string fmt,
    bool total_hist,
    bool virus_info,
    bool virus_hist,
    bool tool_info,
    bool tool_hist,
    bool transmission,
    bool transition,
    bool reproductive,
    bool generation
    )
{

    // Counting number of %
    int n_fmt = 0;
    for (auto & f : fmt)
        if (f == '%')
            n_fmt++;

    if (n_fmt != 1)
        throw std::logic_error("The -fmt- argument must have only one \"%\" symbol.");

    // Listting things to save
    std::vector< bool > what_to_save = {
        virus_info,
        virus_hist,
        tool_info,
        tool_hist,
        total_hist,
        transmission,
        transition,
        reproductive,
        generation
    };

    std::function<void(size_t,Model<TSeq>*)> saver = [fmt,what_to_save](
        size_t niter, Model<TSeq> * m
    ) -> void {

        std::string virus_info = "";
        std::string virus_hist = "";
        std::string tool_info = "";
        std::string tool_hist = "";
        std::string total_hist = "";
        std::string transmission = "";
        std::string transition = "";
        std::string reproductive = "";
        std::string generation = "";

        char buff[1024u];
        if (what_to_save[0u])
        {
            virus_info = fmt + std::string("_virus_info.csv");
            snprintf(buff, sizeof(buff), virus_info.c_str(), niter);
            virus_info = buff;
        } 
        if (what_to_save[1u])
        {
            virus_hist = fmt + std::string("_virus_hist.csv");
            snprintf(buff, sizeof(buff), virus_hist.c_str(), niter);
            virus_hist = buff;
        } 
        if (what_to_save[2u])
        {
            tool_info = fmt + std::string("_tool_info.csv");
            snprintf(buff, sizeof(buff), tool_info.c_str(), niter);
            tool_info = buff;
        } 
        if (what_to_save[3u])
        {
            tool_hist = fmt + std::string("_tool_hist.csv");
            snprintf(buff, sizeof(buff), tool_hist.c_str(), niter);
            tool_hist = buff;
        } 
        if (what_to_save[4u])
        {
            total_hist = fmt + std::string("_total_hist.csv");
            snprintf(buff, sizeof(buff), total_hist.c_str(), niter);
            total_hist = buff;
        } 
        if (what_to_save[5u])
        {
            transmission = fmt + std::string("_transmission.csv");
            snprintf(buff, sizeof(buff), transmission.c_str(), niter);
            transmission = buff;
        } 
        if (what_to_save[6u])
        {
            transition = fmt + std::string("_transition.csv");
            snprintf(buff, sizeof(buff), transition.c_str(), niter);
            transition = buff;
        } 
        if (what_to_save[7u])
        {

            reproductive = fmt + std::string("_reproductive.csv");
            snprintf(buff, sizeof(buff), reproductive.c_str(), niter);
            reproductive = buff;

        }
        if (what_to_save[8u])
        {

            generation = fmt + std::string("_generation.csv");
            snprintf(buff, sizeof(buff), generation.c_str(), niter);
            generation = buff;

        }
        
    
        m->write_data(
            virus_info,
            virus_hist,
            tool_info,
            tool_hist,
            total_hist,
            transmission,
            transition,
            reproductive,
            generation
        );

    };

    return saver;
}


template<typename TSeq>
inline void Model<TSeq>::events_add(
    Agent<TSeq> * agent_,
    VirusPtr<TSeq> virus_,
    ToolPtr<TSeq> tool_,
    Entity<TSeq> * entity_,
    epiworld_fast_int new_state_,
    epiworld_fast_int queue_,
    EventFun<TSeq> call_,
    int idx_agent_,
    int idx_object_
) {

    ++nactions;

    #ifdef EPI_DEBUG
    if (nactions == 0)
        throw std::logic_error("Events cannot be zero!!");
    #endif

    if (nactions > events.size())
    {

        events.emplace_back(
            Event<TSeq>(
                agent_, virus_, tool_, entity_, new_state_, queue_, call_,
                idx_agent_, idx_object_
            ));

    }
    else 
    {

        Event<TSeq> & A = events.at(nactions - 1u);

        A.agent      = std::move(agent_);
        A.virus      = std::move(virus_);
        A.tool       = std::move(tool_);
        A.entity     = std::move(entity_);
        A.new_state  = std::move(new_state_);
        A.queue      = std::move(queue_);
        A.call       = std::move(call_);
        A.idx_agent  = std::move(idx_agent_);
        A.idx_object = std::move(idx_object_);

    }

    return;

}

template<typename TSeq>
inline void Model<TSeq>::events_run()
{
    // Making the call
    size_t nevents_tmp = 0;
    while (nevents_tmp < nactions)
    {

        Event<TSeq> & a = events[nevents_tmp++];
        Agent<TSeq> * p  = a.agent;

        #ifdef EPI_DEBUG
        if (a.new_state >= static_cast<epiworld_fast_int>(nstates))
        {
            throw std::range_error(
                "The proposed state " + std::to_string(a.new_state) + " is out of range. " +
                "The model currently has " + std::to_string(nstates - 1) + " states.");

        }
        else if ((a.new_state != -99) && (a.new_state < 0))
        {
            throw std::range_error(
                "The proposed state " + std::to_string(a.new_state) + " is out of range. " +
                "The state cannot be negative.");
        }
        #endif

        // Undoing the change in the transition matrix
        if (
            (a.new_state != -99) &&
            (p->state_last_changed == today()) &&
            (static_cast<int>(p->state) != a.new_state)
        )
        {
            // Undoing state change in the transition matrix
            // The previous state is already recorded
            db.update_state(p->state_prev, p->state, true);

        } else if (p->state_last_changed != today()) 
            p->state_prev = p->state; // Recording the previous state

        // Applying function after the fact. This way, if there were
        // updates, they can be recorded properly, before losing the information
        if (a.call)
        {
            a.call(a, this);
        }

        if (a.new_state != -99)
            p->state = a.new_state;

        // Registering that the last change was today
        p->state_last_changed = today();
        

        #ifdef EPI_DEBUG
        if (static_cast<int>(p->state) >= static_cast<int>(nstates))
                throw std::range_error(
                    "The new state " + std::to_string(p->state) + " is out of range. " +
                    "The model currently has " + std::to_string(nstates - 1) + " states.");
        #endif

        // Updating queue
        if (use_queuing && a.queue != -99)
        {

            if (a.queue == Queue<TSeq>::Everyone)
                queue += p;
            else if (a.queue == -Queue<TSeq>::Everyone)
                queue -= p;
            else if (a.queue == Queue<TSeq>::OnlySelf)
                queue[p->get_id()]++;
            else if (a.queue == -Queue<TSeq>::OnlySelf)
                queue[p->get_id()]--;
            else if (a.queue != Queue<TSeq>::NoOne)
                throw std::logic_error(
                    "The proposed queue change is not valid. Queue values can be {-2, -1, 0, 1, 2}."
                    );
                    
        }

    }

    // Go back to square 1
    nactions = 0u;

    return;
    
}

/**
 * @name Default function for combining susceptibility_reduction levels
 * 
 * @tparam TSeq 
 * @param pt 
 * @return epiworld_double 
 */
///@{
template<typename TSeq>
inline epiworld_double susceptibility_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq> * m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_susceptibility_reduction(v, m));

    return 1.0 - total;
    
}

template<typename TSeq>
inline epiworld_double transmission_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_transmission_reduction(v, m));

    return (1.0 - total);
    
}

template<typename TSeq>
inline epiworld_double recovery_enhancer_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
)
{
    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
        total *= (1.0 - tool->get_recovery_enhancer(v, m));

    return 1.0 - total;
    
}

template<typename TSeq>
inline epiworld_double death_reduction_mixer_default(
    Agent<TSeq>* p,
    VirusPtr<TSeq> v,
    Model<TSeq>* m
) {

    epiworld_double total = 1.0;
    for (auto & tool : p->get_tools())
    {
        total *= (1.0 - tool->get_death_reduction(v, m));
    } 

    return 1.0 - total;
    
}
///@}

template<typename TSeq>
inline Model<TSeq> * Model<TSeq>::clone_ptr()
{
    // Everything is copied
    Model<TSeq> * ptr = new Model<TSeq>(*dynamic_cast<const Model<TSeq>*>(this));

    #ifdef EPI_DEBUG
    if (*this != *ptr)
        throw std::logic_error("Model::clone_ptr The copies of the model don't match.");
    #endif

    return ptr;
}

template<typename TSeq>
inline Model<TSeq>::Model()
{
    db.model = this;
    db.user_data = this;
    if (use_queuing)
        queue.model = this;
}

template<typename TSeq>
inline Model<TSeq>::Model(const Model<TSeq> & model) :
    name(model.name),
    db(model.db),
    population(model.population),
    population_backup(model.population_backup),
    directed(model.directed),
    viruses(model.viruses),
    tools(model.tools),
    entities(model.entities),
    entities_backup(model.entities_backup),
    rewire_fun(model.rewire_fun),
    rewire_prop(model.rewire_prop),
    parameters(model.parameters),
    ndays(model.ndays),
    pb(model.pb),
    state_fun(model.state_fun),
    states_labels(model.states_labels),
    initial_states_fun(model.initial_states_fun),
    nstates(model.nstates),
    verbose(model.verbose),
    current_date(model.current_date),
    globalevents(model.globalevents),
    queue(model.queue),
    use_queuing(model.use_queuing),
    array_double_tmp(model.array_double_tmp.size()),
    array_virus_tmp(model.array_virus_tmp.size())
{


    // Removing old neighbors
    for (auto & p : population)
        p.model = this;

    // Pointing to the right place. This needs
    // to be done afterwards since the state zero is set as a function
    // of the population.
    db.model = this;
    db.user_data.model = this;

    if (use_queuing)
        queue.model = this;

    agents_data = model.agents_data;
    agents_data_ncols = model.agents_data_ncols;

}

template<typename TSeq>
inline Model<TSeq>::Model(Model<TSeq> & model) :
    Model(dynamic_cast< const Model<TSeq> & >(model)) {}

template<typename TSeq>
inline Model<TSeq>::Model(Model<TSeq> && model) :
    name(std::move(model.name)),
    db(std::move(model.db)),
    population(std::move(model.population)),
    agents_data(std::move(model.agents_data)),
    agents_data_ncols(std::move(model.agents_data_ncols)),
    directed(std::move(model.directed)),
    // Virus
    viruses(std::move(model.viruses)),
    // Tools
    tools(std::move(model.tools)),
    // Entities
    entities(std::move(model.entities)),
    entities_backup(std::move(model.entities_backup)),
    // Pseudo-RNG
    engine(std::move(model.engine)),
    runifd(std::move(model.runifd)),
    rnormd(std::move(model.rnormd)),
    rgammad(std::move(model.rgammad)),
    rlognormald(std::move(model.rlognormald)),
    rexpd(std::move(model.rexpd)),
    // Rewiring
    rewire_fun(std::move(model.rewire_fun)),
    rewire_prop(std::move(model.rewire_prop)),
    parameters(std::move(model.parameters)),
    // Others
    ndays(model.ndays),
    pb(std::move(model.pb)),
    state_fun(std::move(model.state_fun)),
    states_labels(std::move(model.states_labels)),
    initial_states_fun(std::move(model.initial_states_fun)),
    nstates(model.nstates),
    verbose(model.verbose),
    current_date(std::move(model.current_date)),
    globalevents(std::move(model.globalevents)),
    queue(std::move(model.queue)),
    use_queuing(model.use_queuing),
    array_double_tmp(model.array_double_tmp.size()),
    array_virus_tmp(model.array_virus_tmp.size())
{

    db.model = this;
    db.user_data.model = this;

    if (use_queuing)
        queue.model = this;

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::operator=(const Model<TSeq> & m)
{
    name = m.name;

    population        = m.population;
    population_backup = m.population_backup;

    for (auto & p : population)
        p.model = this;

    db = m.db;
    db.model = this;
    db.user_data.model = this;

    directed = m.directed;
    
    viruses                        = m.viruses;

    tools                         = m.tools;
    
    entities        = m.entities;
    entities_backup = m.entities_backup;
    
    rewire_fun  = m.rewire_fun;
    rewire_prop = m.rewire_prop;

    parameters = m.parameters;
    ndays      = m.ndays;
    pb         = m.pb;

    state_fun    = m.state_fun;
    states_labels = m.states_labels;
    initial_states_fun = m.initial_states_fun;
    nstates       = m.nstates;

    verbose     = m.verbose;

    current_date = m.current_date;

    globalevents = m.globalevents;

    queue       = m.queue;
    use_queuing = m.use_queuing;

    // Making sure population is passed correctly
    // Pointing to the right place
    db.model = this;
    db.user_data.model = this;

    agents_data            = m.agents_data;
    agents_data_ncols = m.agents_data_ncols;

    // Figure out the queuing
    if (use_queuing)
        queue.model = this;

    array_double_tmp.resize(static_cast<size_t>(1024u), 0.0);
    array_virus_tmp.resize(1024u);

    return *this;

}

template<typename TSeq>
inline DataBase<TSeq> & Model<TSeq>::get_db()
{
    return db;
}

template<typename TSeq>
inline const DataBase<TSeq> & Model<TSeq>::get_db() const
{
    return db;
}


template<typename TSeq>
inline std::vector<Agent<TSeq>> & Model<TSeq>::get_agents()
{
    return population;
}

template<typename TSeq>
inline Agent<TSeq> & Model<TSeq>::get_agent(size_t i)
{
    return population[i];
}

template<typename TSeq>
inline std::vector< epiworld_fast_uint > Model<TSeq>::get_agents_states() const
{
    std::vector< epiworld_fast_uint > states(population.size());
    for (size_t i = 0u; i < population.size(); ++i)
        states[i] = population[i].get_state();

    return states;
}

template<typename TSeq>
inline std::vector< Viruses_const<TSeq> > Model<TSeq>::get_agents_viruses() const
{

    std::vector< Viruses_const<TSeq> > viruses(population.size());
    for (size_t i = 0u; i < population.size(); ++i)
        viruses[i] = population[i].get_virus();

    return viruses;

}

// Same as before, but the non const version
template<typename TSeq>
inline std::vector< Viruses<TSeq> > Model<TSeq>::get_agents_viruses()
{

    std::vector< Viruses<TSeq> > viruses(population.size());
    for (size_t i = 0u; i < population.size(); ++i)
        viruses[i] = population[i].get_virus();

    return viruses;

}

template<typename TSeq>
inline std::vector<Entity<TSeq>> & Model<TSeq>::get_entities()
{
    return entities;
}

template<typename TSeq>
inline Entity<TSeq> & Model<TSeq>::get_entity(size_t i, int * entity_pos)
{
    
    for (size_t j = 0u; j < entities.size(); ++j)
        if (entities[j].get_id() == static_cast<int>(i))
        {

            if (entity_pos)
                *entity_pos = j;

            return entities[j];

        }

    throw std::range_error("The entity with id " + std::to_string(i) + " was not found.");

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::agents_smallworld(
    epiworld_fast_uint n,
    epiworld_fast_uint k,
    bool d,
    epiworld_double p
)
{
    agents_from_adjlist(
        rgraph_smallworld(n, k, p, d, *this)
    );

    return *this;
}

template<typename TSeq>
inline void Model<TSeq>::agents_empty_graph(
    epiworld_fast_uint n
) 
{

    // Resizing the people
    population.clear();
    population.resize(n, Agent<TSeq>());

    // Filling the model and ids
    size_t i = 0u;
    for (auto & p : population)
    {
        p.id = i++;
        p.model = this;
    }
    

}

template<typename TSeq>
inline void Model<TSeq>::set_rand_gamma(epiworld_double alpha, epiworld_double beta)
{
    rgammad = std::gamma_distribution<>(alpha,beta);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_norm(epiworld_double mean, epiworld_double sd)
{ 
    rnormd  = std::normal_distribution<>(mean, sd);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_unif(epiworld_double a, epiworld_double b)
{ 
    runifd  = std::uniform_real_distribution<>(a, b);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_lognormal(epiworld_double mean, epiworld_double shape)
{ 
    rlognormald  = std::lognormal_distribution<>(mean, shape);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_exp(epiworld_double lambda)
{ 
    rexpd  = std::exponential_distribution<>(lambda);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_binom(int n, epiworld_double p)
{ 
    rbinomd  = std::binomial_distribution<>(n, p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_nbinom(int n, epiworld_double p)
{ 
    rnbinomd  = std::negative_binomial_distribution<>(n, p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_geom(epiworld_double p)
{ 
    rgeomd  = std::geometric_distribution<>(p);
}

template<typename TSeq>
inline void Model<TSeq>::set_rand_poiss(epiworld_double lambda)
{ 
    rpoissd  = std::poisson_distribution<>(lambda);
}

template<typename TSeq>
inline epiworld_double & Model<TSeq>::operator()(std::string pname) {

    if (parameters.find(pname) == parameters.end())
        throw std::range_error("The parameter '"+ pname + "' is not in the model.");

    return parameters[pname];

}

template<typename TSeq>
inline size_t Model<TSeq>::size() const {
    return population.size();
}

template<typename TSeq>
inline void Model<TSeq>::dist_virus()
{

    for (auto & v: viruses)
    {

        v->distribute(this);

        // Apply the events
        events_run();
    }

}

template<typename TSeq>
inline void Model<TSeq>::dist_tools()
{

    for (auto & tool: tools)
    {

        tool->distribute(this);

        // Apply the events
        events_run();

    }

}

template<typename TSeq>
inline void Model<TSeq>::dist_entities()
{

    for (auto & entity: entities)
    {

        entity.distribute(this);

        // Apply the events
        events_run();

    }

}

template<typename TSeq>
inline void Model<TSeq>::chrono_start() {
    time_start = std::chrono::steady_clock::now();
}

template<typename TSeq>
inline void Model<TSeq>::chrono_end() {
    time_end = std::chrono::steady_clock::now();
    time_elapsed += (time_end - time_start);
    n_replicates++;
}

template<typename TSeq>
inline void Model<TSeq>::set_backup()
{

    if (population_backup.size() == 0u)
        population_backup = std::vector< Agent<TSeq> >(population);

    if (entities_backup.size() == 0u)
        entities_backup = std::vector< Entity<TSeq> >(entities);

}

template<typename TSeq>
inline std::shared_ptr< std::mt19937 > & Model<TSeq>::get_rand_endgine()
{
    return engine;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif() {
    // CHECK_INIT()
    return runifd(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::runif(epiworld_double a, epiworld_double b) {
    // CHECK_INIT()
    return runifd(*engine) * (b - a) + a;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm() {
    // CHECK_INIT()
    return rnormd(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rnorm(epiworld_double mean, epiworld_double sd) {
    // CHECK_INIT()
    return rnormd(*engine) * sd + mean;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma() {
    return rgammad(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rgamma(epiworld_double alpha, epiworld_double beta) {

    return rgammad(
        *engine,
        std::gamma_distribution<>::param_type(alpha, beta)
    );
    
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp() {
    return rexpd(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rexp(epiworld_double lambda) {

    return rexpd(
        *engine,
        std::exponential_distribution<>::param_type(lambda)
    );

}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal() {
    return rlognormald(*engine);
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::rlognormal(epiworld_double mean, epiworld_double shape) {
    
    return rlognormald(
        *engine,
        std::lognormal_distribution<>::param_type(mean, shape)
    );
}

template<typename TSeq>
inline int Model<TSeq>::rbinom() {
    return rbinomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rbinom(int n, epiworld_double p) {

    if (n == 0 || p == 0.0)
        return 0;

    return rbinomd(
        *engine,
        std::binomial_distribution<>::param_type(n, p)
    );

}

template<typename TSeq>
inline int Model<TSeq>::rnbinom() {
    return rnbinomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rnbinom(int n, epiworld_double p) {

    return rnbinomd(
        *engine,
        std::negative_binomial_distribution<>::param_type(n, p)
    );
}

template<typename TSeq>
inline int Model<TSeq>::rgeom() {
    return rgeomd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rgeom(epiworld_double p) {

    return rgeomd(
        *engine,
        std::geometric_distribution<>::param_type(p)
    );

}

template<typename TSeq>
inline int Model<TSeq>::rpoiss() {
    return rpoissd(*engine);
}

template<typename TSeq>
inline int Model<TSeq>::rpoiss(epiworld_double lambda) {
    
    return rpoissd(
        *engine,
        std::poisson_distribution<>::param_type(lambda)
    );

}

template<typename TSeq>
inline void Model<TSeq>::seed(size_t s) {
    this->engine->seed(s);
}

template<typename TSeq>
inline void Model<TSeq>::add_virus(
    Virus<TSeq> & v
    )
{

    // Checking the state
    epiworld_fast_int init_, post_, rm_;
    v.get_state(&init_, &post_, &rm_);

    if (init_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -init- state."
            );
    else if (post_ == -99)
        throw std::logic_error(
            "The virus \"" + v.get_name() + "\" has no -post- state."
            );
    
    // Recording the variant
    db.record_virus(v);

    // Adding new virus
    viruses.push_back(std::make_shared< Virus<TSeq> >(v));

}

template<typename TSeq>
inline void Model<TSeq>::add_tool(Tool<TSeq> & t)
{

    
    db.record_tool(t);

    // Adding the tool to the model (and database.)
    tools.push_back(std::make_shared< Tool<TSeq> >(t));

}

template<typename TSeq>
inline void Model<TSeq>::add_entity(Entity<TSeq> e)
{

    e.id = entities.size();
    entities.push_back(e);

}

template<typename TSeq>
inline void Model<TSeq>::rm_entity(size_t entity_id)
{

    int entity_pos = 0;
    auto & entity = this->get_entity(entity_id, &entity_pos);

    // First, resetting the entity
    entity.reset();

    // How should
    if (entity_pos != (static_cast<int>(entities.size()) - 1))
        std::swap(entities[entity_pos], entities[entities.size() - 1]);

    entities.pop_back();
}

template<typename TSeq>
inline void Model<TSeq>::rm_virus(size_t virus_pos)
{

    if (viruses.size() <= virus_pos)
        throw std::range_error(
            std::string("The specified virus (") +
            std::to_string(virus_pos) +
            std::string(") is out of range. ") +
            std::string("There are only ") +
            std::to_string(viruses.size()) +
            std::string(" viruses.")
            );

    // Flipping with the last one
    std::swap(viruses[virus_pos], viruses[viruses.size() - 1]);
    viruses.pop_back();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::rm_tool(size_t tool_pos)
{

    if (tools.size() <= tool_pos)
        throw std::range_error(
            std::string("The specified tool (") +
            std::to_string(tool_pos) +
            std::string(") is out of range. ") +
            std::string("There are only ") +
            std::to_string(tools.size()) +
            std::string(" tools.")
            );

    // Flipping with the last one
    std::swap(tools[tool_pos], tools[tools.size() - 1]);
    
    /* There's an error on windows:
    https://github.com/UofUEpiBio/epiworldR/actions/runs/4801482395/jobs/8543744180#step:6:84

    More clear here:
    https://stackoverflow.com/questions/58660207/why-doesnt-stdswap-work-on-vectorbool-elements-under-clang-win
    */

    tools.pop_back();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::load_agents_entities_ties(
    std::string fn,
    int skip
    )
{

    int i,j;
    std::ifstream filei(fn);

    if (!filei)
        throw std::logic_error("The file " + fn + " was not found.");

    int linenum = 0;
    std::vector< std::vector< epiworld_fast_uint > > target_(entities.size());

    target_.reserve(1e5);

    while (!filei.eof())
    {

        if (linenum++ < skip)
            continue;

        filei >> i >> j;

        // Looking for exceptions
        if (filei.bad())
            throw std::logic_error(
                "I/O error while reading the file " +
                fn
            );

        if (filei.fail())
            break;

        if (i >= static_cast<int>(this->size()))
            throw std::range_error(
                "The agent["+std::to_string(linenum)+"] = " + std::to_string(i) +
                " is above the max id " + std::to_string(this->size() - 1)
                );

        if (j >= static_cast<int>(this->entities.size()))
            throw std::range_error(
                "The entity["+std::to_string(linenum)+"] = " + std::to_string(j) +
                " is above the max id " + std::to_string(this->entities.size() - 1)
                );

        target_[j].push_back(i);

        population[i].add_entity(entities[j], nullptr);

    }

    return;

}

template<typename TSeq>
inline void Model<TSeq>::load_agents_entities_ties(
    const std::vector< int > & agents_ids,
    const std::vector< int > & entities_ids
) {

    // Checking the size
    if (agents_ids.size() != entities_ids.size())
        throw std::length_error(
            std::string("The size of agents_ids (") +
            std::to_string(agents_ids.size()) +
            std::string(") and entities_ids (") +
            std::to_string(entities_ids.size()) +
            std::string(") must be the same.")
            );

    return this->load_agents_entities_ties(
        agents_ids.data(),
        entities_ids.data(),
        agents_ids.size()
    );

}

template<typename TSeq>
inline void Model<TSeq>::load_agents_entities_ties(
    const int * agents_ids,
    const int * entities_ids,
    size_t n
) {

    auto get_agent = [agents_ids](int i) -> int {
        return *(agents_ids + i);
        };

    auto get_entity = [entities_ids](int i) -> int {
        return *(entities_ids + i);
        };

    for (size_t i = 0u; i < n; ++i)
    {

        if (get_agent(i) < 0)
            throw std::length_error(
                std::string("agents_ids[") +
                std::to_string(i) +
                std::string("] = ") +
                std::to_string(get_agent(i)) +
                std::string(" is negative.")
                );

        if (get_entity(i) < 0)
            throw std::length_error(
                std::string("entities_ids[") +
                std::to_string(i) +
                std::string("] = ") +
                std::to_string(get_entity(i)) +
                std::string(" is negative.")
                );

        int pop_size = static_cast<int>(this->population.size());
        if (get_agent(i) >= pop_size)
            throw std::length_error(
                std::string("agents_ids[") +
                std::to_string(i) +
                std::string("] = ") +
                std::to_string(get_agent(i)) +
                std::string(" is out of range (population size: ") +
                std::to_string(pop_size) +
                std::string(").")
                );

        int ent_size = static_cast<int>(this->entities.size());
        if (get_entity(i) >= ent_size)
            throw std::length_error(
                std::string("entities_ids[") +
                std::to_string(i) +
                std::string("] = ") +
                std::to_string(get_entity(i)) +
                std::string(" is out of range (entities size: ") +
                std::to_string(ent_size) +
                std::string(").")
                );

        // Adding the entity to the agent
        this->population[get_agent(i)].add_entity(
            this->entities[get_entity(i)],
            nullptr /* Immediately add it to the agent */
        );

    }

    return;


}

template<typename TSeq>
inline void Model<TSeq>::agents_from_adjlist(
    std::string fn,
    int size,
    int skip,
    bool directed
    ) {

    AdjList al;
    al.read_edgelist(fn, size, skip, directed);
    this->agents_from_adjlist(al);

}

template<typename TSeq>
inline void Model<TSeq>::agents_from_edgelist(
    const std::vector< int > & source,
    const std::vector< int > & target,
    int size,
    bool directed
) {

    AdjList al(source, target, size, directed);
    agents_from_adjlist(al);

}

template<typename TSeq>
inline void Model<TSeq>::agents_from_adjlist(AdjList al) {

    // Resizing the people
    agents_empty_graph(al.vcount());
    
    const auto & tmpdat = al.get_dat();
    
    for (size_t i = 0u; i < tmpdat.size(); ++i)
    {

        // population[i].id    = i;
        population[i].model = this;

        for (const auto & link: tmpdat[i])
        {

            population[i].add_neighbor(
                population[link.first],
                true, true
                );

        }

    }

    #ifdef EPI_DEBUG
    for (auto & p: population)
    {
        if (p.id >= static_cast<int>(al.vcount()))
            throw std::logic_error(
                "Agent's id cannot be negative above or equal to the number of agents!");
    }
    #endif

}

template<typename TSeq>
inline bool Model<TSeq>::is_directed() const
{
    if (population.size() == 0u)
        throw std::logic_error("The population hasn't been initialized.");

    return directed;
}

template<typename TSeq>
inline int Model<TSeq>::today() const {

    if (ndays == 0)
      return 0;

    return this->current_date;
}

template<typename TSeq>
inline void Model<TSeq>::next() {

    #ifdef EPI_DEBUG
    // Checking all the agents have proper states
    for (auto & p : population)
    {
        if ((p.state >= nstates) || (p.state < 0))
        {
            throw std::range_error(
                "The agent " + std::to_string(p.id) +
                " has state " + std::to_string(p.state) +
                " which is above the maximum state of " +
                std::to_string(nstates - 1) + "."
            );
        }

        if ((p.state_prev >= nstates) || (p.state_prev < 0))
        {
            throw std::range_error(
                "The agent " + std::to_string(p.id) +
                " has previous state " + std::to_string(p.state_prev) +
                " which is above the maximum state of " +
                std::to_string(nstates - 1) + "."
            );
        }
        
    }

    #endif

    db.record();
    ++this->current_date;
    
    // Advancing the progress bar
    if ((this->current_date >= 1) && verbose)
        pb.next();

    return ;
}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::run(
    epiworld_fast_uint ndays,
    int seed
) 
{

    if (size() == 0u)
        throw std::logic_error("There are no agents in this model!");

    if (nstates == 0u)
        throw std::logic_error(
            std::string("No states registered in this model. ") +
            std::string("At least one state should be included. See the function -Model::add_state()-")
            );

    // Setting up the number of steps
    this->ndays = ndays;

    if (seed >= 0)
        engine->seed(seed);

    array_double_tmp.resize(std::max(
        size(),
        static_cast<size_t>(1024)
    ));


    array_virus_tmp.resize(1024);

    // Checking whether the proposed state in/out/removed
    // are valid
    epiworld_fast_int _init, _end, _removed;
    int nstate_int = static_cast<int>(nstates);
    for (auto & v : viruses)
    {
        v->get_state(&_init, &_end, &_removed);
        
        // Negative unspecified state
        if (((_init != -99) && (_init < 0)) || (_init >= nstate_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstates - 1));

        // Negative unspecified state
        if (((_end != -99) && (_end < 0)) || (_end >= nstate_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstates - 1));

        if (((_removed != -99) && (_removed < 0)) || (_removed >= nstate_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstates - 1));

    }

    for (auto & t : tools)
    {
        t->get_state(&_init, &_end);
        
        // Negative unspecified state
        if (((_init != -99) && (_init < 0)) || (_init >= nstate_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstates - 1));

        // Negative unspecified state
        if (((_end != -99) && (_end < 0)) || (_end >= nstate_int))
            throw std::range_error("States must be between 0 and " +
                std::to_string(nstates - 1));

    }

    // Starting first infection and tools
    reset();

    // Initializing the simulation
    chrono_start();
    EPIWORLD_RUN((*this))
    {

        #ifdef EPI_DEBUG
        db.n_transmissions_potential = 0;
        db.n_transmissions_today = 0;
        #endif

        // We can execute these components in whatever order the
        // user needs.
        this->update_state();
    
        // We start with the Global events
        this->run_globalevents();

        // In this case we are applying degree sequence rewiring
        // to change the network just a bit.
        this->rewire();

        // This locks all the changes
        this->next();

        // Mutation must happen at the very end of all
        this->mutate_virus();

    }

    // The last reaches the end...
    this->current_date--;

    chrono_end();

    return *this;

}

template<typename TSeq>
inline void Model<TSeq>::run_multiple(
    epiworld_fast_uint ndays,
    epiworld_fast_uint nexperiments,
    int seed_,
    std::function<void(size_t,Model<TSeq>*)> fun,
    bool reset,
    bool verbose,
    #ifdef _OPENMP
    int nthreads
    #else
    int
    #endif
)
{

    if (seed_ >= 0)
        this->seed(seed_);

    // Seeds will be reproducible by default
    std::vector< int > seeds_n(nexperiments);
    for (auto & s : seeds_n)
    {
        s = static_cast<int>(
            std::floor(
                runif() * static_cast<double>(std::numeric_limits<int>::max())
                )
        );
    }
    // #endif

    EPI_DEBUG_NOTIFY_ACTIVE()

    bool old_verb = this->verbose;
    verbose_off();

    // Setting up backup
    if (reset)
        set_backup();

    #ifdef _OPENMP

    omp_set_num_threads(nthreads);

    // Not more than the number of experiments
    nthreads =
        static_cast<size_t>(nthreads) > nexperiments ? nexperiments : nthreads;

    // Generating copies of the model
    std::vector< Model<TSeq> * > these(
        std::max(nthreads - 1, 0)
    );

    #pragma omp parallel for shared(these, nthreads)
    for (size_t i = 0u; i < static_cast<size_t>(nthreads); ++i)
    {

        if (i == 0)
            continue;

        these[i - 1] = clone_ptr();

    }
        

    // Figuring out how many replicates - distribute remainder evenly
    std::vector< size_t > nreplicates(nthreads, 0);
    std::vector< size_t > nreplicates_csum(nthreads, 0);
    
    size_t base_replicates = nexperiments / nthreads;
    size_t remainder = nexperiments % nthreads;
    
    size_t sums = 0u;
    for (int i = 0; i < nthreads; ++i)
    {
        // Distribute remainder to first 'remainder' threads
        nreplicates[i] = base_replicates + (static_cast<size_t>(i) < remainder ? 1 : 0);
        
        // This takes the cumsum
        nreplicates_csum[i] = sums;
        sums += nreplicates[i];
    }

    Progress pb_multiple(
        nreplicates[0u],
        EPIWORLD_PROGRESS_BAR_WIDTH
        );

    if (verbose)
    {

        printf_epiworld(
            "Starting multiple runs (%i) using %i thread(s)\n", 
            static_cast<int>(nexperiments),
            static_cast<int>(nthreads)
        );

        pb_multiple.start();

    }

    #ifdef EPI_DEBUG
    // Checking the initial state of all the models. Throw an
    // exception if they are not the same.
    for (size_t i = 1; i < static_cast<size_t>(std::max(nthreads - 1, 0)); ++i)
    {

        if (db != these[i]->db)
        {
            throw std::runtime_error(
                "The initial state of the models is not the same"
            );
        }
    }
    #endif
    
    #pragma omp parallel shared(these, nreplicates, nreplicates_csum, seeds_n) \
        firstprivate(nexperiments, nthreads, fun, reset, verbose, pb_multiple, ndays) \
        default(shared)
    {

        auto iam = omp_get_thread_num();

        for (size_t n = 0u; n < nreplicates[iam]; ++n)
        {
            size_t sim_id = nreplicates_csum[iam] + n;
            if (iam == 0)
            {

                // Checking if the user interrupted the simulation
                EPI_CHECK_USER_INTERRUPT(n);

                // Initializing the seed
                run(ndays, seeds_n[sim_id]);

                if (fun)
                    fun(n, this);

                // Only the first one prints
                if (verbose)
                    pb_multiple.next();                

            } else {

                // Initializing the seed
                these[iam - 1]->run(ndays, seeds_n[sim_id]);

                if (fun)
                    fun(sim_id, these[iam - 1]);

            }

            

        }
        
    }

    // Adjusting the number of replicates
    n_replicates += (nexperiments - nreplicates[0u]);

    #pragma omp parallel for shared(these)
    for (int i = 1; i < nthreads; ++i)
    {
        delete these[i - 1];
    }

    #else
    // if (reset)
    //     set_backup();

    Progress pb_multiple(
        nexperiments,
        EPIWORLD_PROGRESS_BAR_WIDTH
        )
        ;
    if (verbose)
    {

        printf_epiworld(
            "Starting multiple runs (%i)\n", 
            static_cast<int>(nexperiments)
        );

        pb_multiple.start();

    }

    for (size_t n = 0u; n < nexperiments; ++n)
    {

        // Checking if the user interrupted the simulation
        EPI_CHECK_USER_INTERRUPT(n);

        run(ndays, seeds_n[n]);

        if (fun)
            fun(n, this);

        if (verbose)
            pb_multiple.next();
    
    }
    #endif

    if (verbose)
        pb_multiple.end();

    if (old_verb)
        verbose_on();

    return;

}

template<typename TSeq>
inline void Model<TSeq>::update_state() {

    // Next state
    if (use_queuing)
    {
        int i = -1;
        for (auto & p: population)
            if (queue[++i] > 0)
            {
                if (state_fun[p.state])
                    state_fun[p.state](&p, this);
            }

    }
    else
    {

        for (auto & p: population)
            if (state_fun[p.state])
                    state_fun[p.state](&p, this);

    }

    events_run();
    
}

template<typename TSeq>
inline void Model<TSeq>::mutate_virus() {

    // Checking if any virus has mutation
    size_t nmutates = 0u;
    for (const auto & v: viruses)
        if (v->virus_functions->mutation)
            nmutates++;

    if (nmutates == 0u)
        return;

    if (use_queuing)
    {

        int i = -1;
        for (auto & p: population)
        {

            if (queue[++i] == 0)
                continue;

            if (p.virus != nullptr)
                p.virus->mutate(this);

        }

    }
    else 
    {

        for (auto & p: population)
        {

            if (p.virus != nullptr)
                p.virus->mutate(this);

        }

    }
    

}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_viruses() const {
    return db.size();
}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_tools() const {
    return tools.size();
}

template<typename TSeq>
inline epiworld_fast_uint Model<TSeq>::get_ndays() const {
    return ndays;
}

template<typename TSeq>
inline epiworld_fast_uint Model<TSeq>::get_n_replicates() const
{
    return n_replicates;
}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_entities() const {
    return entities.size();
}

template<typename TSeq>
inline void Model<TSeq>::set_ndays(epiworld_fast_uint ndays) {
    this->ndays = ndays;
}

template<typename TSeq>
inline bool Model<TSeq>::get_verbose() const {
    return verbose;
}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::verbose_on() {
    verbose = true;
    return *this;
}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::verbose_off() {
    verbose = false;
    return *this;
}

template<typename TSeq>
inline void Model<TSeq>::set_rewire_fun(
    std::function<void(std::vector<Agent<TSeq>>*,Model<TSeq>*,epiworld_double)> fun
    ) {
    rewire_fun = fun;
}

template<typename TSeq>
inline void Model<TSeq>::set_rewire_prop(epiworld_double prop)
{

    if (prop < 0.0)
        throw std::range_error("Proportions cannot be negative.");

    if (prop > 1.0)
        throw std::range_error("Proportions cannot be above 1.0.");

    rewire_prop = prop;
}

template<typename TSeq>
inline epiworld_double Model<TSeq>::get_rewire_prop() const {
    return rewire_prop;
}

template<typename TSeq>
inline void Model<TSeq>::rewire() {

    if (rewire_fun)
        rewire_fun(&population, this, rewire_prop);
}


template<typename TSeq>
inline void Model<TSeq>::write_data(
    std::string fn_virus_info,
    std::string fn_virus_hist,
    std::string fn_tool_info,
    std::string fn_tool_hist,
    std::string fn_total_hist,
    std::string fn_transmission,
    std::string fn_transition,
    std::string fn_reproductive_number,
    std::string fn_generation_time
    ) const
{

    db.write_data(
        fn_virus_info, fn_virus_hist,
        fn_tool_info, fn_tool_hist,
        fn_total_hist, fn_transmission, fn_transition,
        fn_reproductive_number, fn_generation_time
        );

}

template<typename TSeq>
inline void Model<TSeq>::write_edgelist(
    std::string fn
    ) const
{

    // Figuring out the writing sequence
    std::vector< const Agent<TSeq> * > wseq(size());
    for (const auto & p: population)
        wseq[p.id] = &p;

    std::ofstream efile(fn, std::ios_base::out);
    efile << "source target\n";
    if (this->is_directed())
    {

        for (const auto & p : wseq)
        {

            if (p->neighbors == nullptr)
                continue;

            for (auto & n : *p->neighbors)
                efile << p->id << " " << n << "\n";
        }

    } else {

        for (const auto & p : wseq)
        {

            if (p->neighbors == nullptr)
                continue;

            for (auto & n : *p->neighbors)
                if (static_cast<int>(p->id) <= static_cast<int>(n))
                    efile << p->id << " " << n << "\n";
        }

    }

}

template<typename TSeq>
inline void Model<TSeq>::write_edgelist(
std::vector< int > & source,
std::vector< int > & target
) const {

    // Figuring out the writing sequence
    std::vector< const Agent<TSeq> * > wseq(size());
    for (const auto & p: population)
        wseq[p.id] = &p;

    if (this->is_directed())
    {

        for (const auto & p : wseq)
        {
            if (p->neighbors == nullptr)
                continue;

            for (auto & n : *p->neighbors)
            {
                source.push_back(static_cast<int>(p->id));
                target.push_back(static_cast<int>(n));
            }
        }

    } else {

        for (const auto & p : wseq)
        {

            if (p->neighbors == nullptr)
                continue;

            for (auto & n : *p->neighbors) {
                if (static_cast<int>(p->id) <= static_cast<int>(n)) {
                    source.push_back(static_cast<int>(p->id));
                    target.push_back(static_cast<int>(n));
                }
            }
        }

    }


}

template<typename TSeq>
inline std::map<std::string,epiworld_double> & Model<TSeq>::params()
{
    return parameters;
}

template<typename TSeq>
inline void Model<TSeq>::reset() {

    // Restablishing people
    pb = Progress(ndays, 80);

    if (population_backup.size())
    {
        population = population_backup;
    
        // Ensuring the population is poiting to the model
        for (auto & p : population)
            p.model = this;

        #ifdef EPI_DEBUG
        for (size_t i = 0; i < population.size(); ++i)
        {

            if (population[i] != (population_backup)[i])
                throw std::logic_error("Model::reset population doesn't match.");

        }
        #endif

    }

    for (auto & p : population)
        p.reset();

    #ifdef EPI_DEBUG
    for (auto & a: population)
    {
        if (a.get_state() != 0u)
            throw std::logic_error("Model::reset population doesn't match."
                "Some agents are not in the baseline state.");
    }
    #endif
        
    if (entities_backup.size())
    {
        entities = entities_backup;

        #ifdef EPI_DEBUG
        for (size_t i = 0; i < entities.size(); ++i)
        {

            if (entities[i] != (entities_backup)[i])
                throw std::logic_error("Model::reset entities don't match.");

        }
        #endif
        
    }

    for (auto & e: entities)
        e.reset();
    
    current_date = 0;

    db.reset();

    // This also clears the queue
    if (use_queuing)
        queue.reset();

    // Re distributing tools and virus
    dist_virus();
    dist_tools();
    dist_entities();

    // Distributing initial state, if specified
    initial_states_fun(this);

    // Recording the original state (at time 0) and advancing
    // to time 1
    next();


}

// Too big to keep here
#include "model-meat-print.hpp"


template<typename TSeq>
inline void Model<TSeq>::add_state(
    std::string lab, 
    UpdateFun<TSeq> fun
)
{

    // Checking it doesn't match
    for (auto & s : states_labels)
        if (s == lab)
            throw std::logic_error("state \"" + s + "\" already registered.");

    states_labels.push_back(lab);
    state_fun.push_back(fun);
    nstates++;

}


template<typename TSeq>
inline const std::vector< std::string > &
Model<TSeq>::get_states() const
{
    return states_labels;
}

template<typename TSeq>
inline size_t Model<TSeq>::get_n_states() const
{
    return nstates;
}

template<typename TSeq>
inline const std::vector< UpdateFun<TSeq> > &
Model<TSeq>::get_state_fun() const
{
    return state_fun;
}

template<typename TSeq>
inline void Model<TSeq>::print_state_codes() const
{

    // Horizontal line
    std::string line = "";
    for (epiworld_fast_uint i = 0u; i < 80u; ++i)
        line += "_";

    printf_epiworld("\n%s\nstates CODES\n\n", line.c_str());

    epiworld_fast_uint nchar = 0u;
    for (auto & p : states_labels)
        if (p.length() > nchar)
            nchar = p.length();
    
    std::string fmt = " %2i = %-" + std::to_string(nchar + 1 + 4) + "s\n";
    for (epiworld_fast_uint i = 0u; i < nstates; ++i)
    {

        printf_epiworld(
            fmt.c_str(),
            i,
            (states_labels[i] + " (S)").c_str()
        );

    }

}



template<typename TSeq>
inline epiworld_double Model<TSeq>::add_param(
    epiworld_double initial_value,
    std::string pname,
    bool overwrite
    ) {

    if (parameters.find(pname) == parameters.end())
        parameters[pname] = initial_value;
    else if (!overwrite)
        throw std::logic_error("The parameter " + pname + " already exists.");
    else
        parameters[pname] = initial_value;
    
    return initial_value;

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::read_params(std::string fn, bool overwrite)
{

    auto params_map = read_yaml<epiworld_double>(fn);

    for (auto & p : params_map)
        add_param(p.second, p.first, overwrite);

    return *this;

}

template<typename TSeq>
inline epiworld_double Model<TSeq>::get_param(std::string pname)
{
    if (parameters.find(pname) == parameters.end())
        throw std::logic_error("The parameter " + pname + " does not exists.");

    return parameters[pname];
}

template<typename TSeq>
inline void Model<TSeq>::set_param(std::string pname, epiworld_double value)
{
    if (parameters.find(pname) == parameters.end())
        throw std::logic_error("The parameter " + pname + " does not exists.");

    parameters[pname] = value;

    return;

}

// // Same as before but using the size_t method
// template<typename TSeq>
// inline void Model<TSeq>::set_param(size_t k, epiworld_double value)
// {
//     if (k >= parameters.size())
//         throw std::logic_error("The parameter index " + std::to_string(k) + " does not exists.");

//     // Access the k-th element of the std::unordered_map parameters


//     *(parameters.begin() + k) = value;

//     return;
// }

template<typename TSeq>
inline epiworld_double Model<TSeq>::par(std::string pname) const
{
    const auto iter = parameters.find(pname);
    if (iter == parameters.end())
        throw std::logic_error("The parameter " + pname + " does not exists.");
    return iter->second;
}

#define DURCAST(tunit,txtunit) {\
        elapsed       = std::chrono::duration_cast<std::chrono:: tunit>(\
            time_end - time_start).count(); \
        elapsed_total = std::chrono::duration_cast<std::chrono:: tunit>(time_elapsed).count(); \
        abbr_unit     = txtunit;}

template<typename TSeq>
inline void Model<TSeq>::get_elapsed(
    std::string unit,
    epiworld_double * last_elapsed,
    epiworld_double * total_elapsed,
    std::string * unit_abbr,
    bool print
) const {

    // Preparing the result
    epiworld_double elapsed, elapsed_total;
    std::string abbr_unit;

    // Figuring out the length
    if (unit == "auto")
    {

        size_t tlength = std::to_string(
            static_cast<int>(floor(time_elapsed.count()))
            ).length();
        
        if (tlength <= 1)
            unit = "nanoseconds";
        else if (tlength <= 3)
            unit = "microseconds";
        else if (tlength <= 6)
            unit = "milliseconds";
        else if (tlength <= 8)
            unit = "seconds";
        else if (tlength <= 9)
            unit = "minutes";
        else 
            unit = "hours";

    }

    if (unit == "nanoseconds")       DURCAST(nanoseconds,"ns")
    else if (unit == "microseconds") DURCAST(microseconds,"\xC2\xB5s")
    else if (unit == "milliseconds") DURCAST(milliseconds,"ms")
    else if (unit == "seconds")      DURCAST(seconds,"s")
    else if (unit == "minutes")      DURCAST(minutes,"m")
    else if (unit == "hours")        DURCAST(hours,"h")
    else
        throw std::range_error("The time unit " + unit + " is not supported.");


    if (last_elapsed != nullptr)
        *last_elapsed = elapsed;
    if (total_elapsed != nullptr)
        *total_elapsed = elapsed_total;
    if (unit_abbr != nullptr)
        *unit_abbr = abbr_unit;

    if (!print)
        return;

    if (n_replicates > 1u)
    {
        printf_epiworld("last run elapsed time : %.2f%s\n",
            elapsed, abbr_unit.c_str());
        printf_epiworld("total elapsed time    : %.2f%s\n",
            elapsed_total, abbr_unit.c_str());
        printf_epiworld("total runs            : %i\n",
            static_cast<int>(n_replicates));
        printf_epiworld("mean run elapsed time : %.2f%s\n",
            elapsed_total/static_cast<epiworld_double>(n_replicates), abbr_unit.c_str());

    } else {
        printf_epiworld("last run elapsed time : %.2f%s.\n", elapsed, abbr_unit.c_str());
    }
}

template<typename TSeq>
inline void Model<TSeq>::set_user_data(std::vector< std::string > names)
{
    db.set_user_data(names);
}

template<typename TSeq>
inline void Model<TSeq>::add_user_data(epiworld_fast_uint j, epiworld_double x)
{
    db.add_user_data(j, x);
}

template<typename TSeq>
inline void Model<TSeq>::add_user_data(std::vector<epiworld_double> x)
{
    db.add_user_data(x);
}

template<typename TSeq>
inline UserData<TSeq> & Model<TSeq>::get_user_data()
{
    return db.get_user_data();
}

template<typename TSeq>
inline void Model<TSeq>::add_globalevent(
    std::function<void(Model<TSeq>*)> fun,
    std::string name,
    int date
)
{

    globalevents.push_back(
        GlobalEvent<TSeq>(
            fun,
            name,
            date
            )
    );

}

template<typename TSeq>
inline void Model<TSeq>::add_globalevent(
    GlobalEvent<TSeq> action
)
{
    globalevents.push_back(action);
}

template<typename TSeq>
GlobalEvent<TSeq> & Model<TSeq>::get_globalevent(
    std::string name
)
{

    for (auto & a : globalevents)
        if (a.name == name)
            return a;

    throw std::logic_error("The global action " + name + " was not found.");

}

template<typename TSeq>
GlobalEvent<TSeq> & Model<TSeq>::get_globalevent(
    size_t index
)
{

    if (index >= globalevents.size())
        throw std::range_error("The index " + std::to_string(index) + " is out of range.");

    return globalevents[index];

}

// Remove implementation
template<typename TSeq>
inline void Model<TSeq>::rm_globalevent(
    std::string name
)
{

    for (auto it = globalevents.begin(); it != globalevents.end(); ++it)
    {
        if (it->get_name() == name)
        {
            globalevents.erase(it);
            return;
        }
    }

    throw std::logic_error("The global action " + name + " was not found.");

}

// Same as above, but the index implementation
template<typename TSeq>
inline void Model<TSeq>::rm_globalevent(
    size_t index
)
{

    if (index >= globalevents.size())
        throw std::range_error("The index " + std::to_string(index) + " is out of range.");

    globalevents.erase(globalevents.begin() + index);

}

template<typename TSeq>
inline void Model<TSeq>::run_globalevents()
{

    for (auto & action: globalevents)
    {
        action(this, today());
        events_run();
    }

}

template<typename TSeq>
inline void Model<TSeq>::queuing_on()
{
    use_queuing = true;
}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::queuing_off()
{
    use_queuing = false;
    return *this;
}

template<typename TSeq>
inline bool Model<TSeq>::is_queuing_on() const
{
    return use_queuing;
}

template<typename TSeq>
inline Queue<TSeq> & Model<TSeq>::get_queue()
{
    return queue;
}

template<typename TSeq>
inline const std::vector< VirusPtr<TSeq> > & Model<TSeq>::get_viruses() const
{
    return viruses;
}

template<typename TSeq>
const std::vector< ToolPtr<TSeq> > & Model<TSeq>::get_tools() const
{
    return tools;
}

template<typename TSeq>
inline Virus<TSeq> & Model<TSeq>::get_virus(size_t id)
{

    if (viruses.size() <= id)
        throw std::length_error("The specified id for the virus is out of range");

    return *viruses[id];

}

template<typename TSeq>
inline Tool<TSeq> & Model<TSeq>::get_tool(size_t id)
{

    if (tools.size() <= id)
        throw std::length_error("The specified id for the tools is out of range");

    return *tools[id];

}


template<typename TSeq>
inline void Model<TSeq>::set_agents_data(double * data_, size_t ncols_)
{
    agents_data = data_;
    agents_data_ncols = ncols_;
}

template<typename TSeq>
inline double * Model<TSeq>::get_agents_data() {
    return this->agents_data;
}

template<typename TSeq>
inline size_t Model<TSeq>::get_agents_data_ncols() const {
    return this->agents_data_ncols;
}


template<typename TSeq>
inline void Model<TSeq>::set_name(std::string name)
{
    this->name = name;
}

template<typename TSeq>
inline std::string Model<TSeq>::get_name() const 
{
    return this->name;
}

#define VECT_MATCH(a, b, c) \
    EPI_DEBUG_FAIL_AT_TRUE(a.size() != b.size(), c) \
    for (size_t __i = 0u; __i < a.size(); ++__i) \
    {\
        EPI_DEBUG_FAIL_AT_TRUE(a[__i] != b[__i], c) \
    }

template<typename TSeq>
inline bool Model<TSeq>::operator==(const Model<TSeq> & other) const
{
    EPI_DEBUG_FAIL_AT_TRUE(name != other.name, "names don't match")
    EPI_DEBUG_FAIL_AT_TRUE(db != other.db, "database don't match")

    VECT_MATCH(population, other.population, "population doesn't match")

    EPI_DEBUG_FAIL_AT_TRUE(
        using_backup != other.using_backup,
        "Model:: using_backup don't match"
        )
    
    if ((population_backup.size() != 0) & (other.population_backup.size() != 0))
    {

        // False is population_backup.size() != other.population_backup.size()
        if (population_backup.size() != other.population_backup.size())
            return false;

        for (size_t i = 0u; i < population_backup.size(); ++i)
        {
            if (population_backup[i] != other.population_backup[i])
                return false;
        }
        
    } else if ((population_backup.size() == 0) & (other.population_backup.size() != 0)) {
        return false;
    } else if ((population_backup.size() != 0) & (other.population_backup.size() == 0))
    {
        return false;
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        agents_data != other.agents_data,
        "Model:: agents_data don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        agents_data_ncols != other.agents_data_ncols,
        "Model:: agents_data_ncols don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        directed != other.directed,
        "Model:: directed don't match"
    )
    
    // Viruses -----------------------------------------------------------------
    EPI_DEBUG_FAIL_AT_TRUE(
        viruses.size() != other.viruses.size(),
        "Model:: viruses.size() don't match"
        )

    for (size_t i = 0u; i < viruses.size(); ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *viruses[i] != *other.viruses[i],
            "Model:: *viruses[i] don't match"
        )
            
    }
    
    // Tools -------------------------------------------------------------------
    EPI_DEBUG_FAIL_AT_TRUE(
        tools.size() != other.tools.size(),
        "Model:: tools.size() don't match"
        )
        
    for (size_t i = 0u; i < tools.size(); ++i)
    {
        EPI_DEBUG_FAIL_AT_TRUE(
            *tools[i] != *other.tools[i],
            "Model:: *tools[i] don't match"
        )
            
    }
    
    VECT_MATCH(
        entities,
        other.entities,
        "entities don't match"
    )

    if ((entities_backup.size() != 0) & (other.entities_backup.size() != 0))
    {
        
        for (size_t i = 0u; i < entities_backup.size(); ++i)
        {

            EPI_DEBUG_FAIL_AT_TRUE(
                entities_backup[i] != other.entities_backup[i],
                "Model:: entities_backup[i] don't match"
            )

        }
        
    } else if ((entities_backup.size() == 0) & (other.entities_backup.size() != 0)) {
        EPI_DEBUG_FAIL_AT_TRUE(true, "entities_backup don't match")
    } else if ((entities_backup.size() != 0) & (other.entities_backup.size() == 0))
    {
        EPI_DEBUG_FAIL_AT_TRUE(true, "entities_backup don't match")
    }

    EPI_DEBUG_FAIL_AT_TRUE(
        rewire_prop != other.rewire_prop,
        "Model:: rewire_prop don't match"
    )
        
    EPI_DEBUG_FAIL_AT_TRUE(
        parameters.size() != other.parameters.size(),
        "Model:: () don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        parameters != other.parameters,
        "Model:: parameters don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        ndays != other.ndays,
        "Model:: ndays don't match"
    )
    
    VECT_MATCH(
        states_labels,
        other.states_labels,
        "state labels don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        nstates != other.nstates,
        "Model:: nstates don't match"
    )
    
    EPI_DEBUG_FAIL_AT_TRUE(
        verbose != other.verbose,
        "Model:: verbose don't match"
    )

    EPI_DEBUG_FAIL_AT_TRUE(
        current_date != other.current_date,
        "Model:: current_date don't match"
    )

    VECT_MATCH(globalevents, other.globalevents, "global action don't match");

    EPI_DEBUG_FAIL_AT_TRUE(
        queue != other.queue,
        "Model:: queue don't match"
    )
    

    EPI_DEBUG_FAIL_AT_TRUE(
        use_queuing != other.use_queuing,
        "Model:: use_queuing don't match"
    )
    
    return true;

}

template<typename TSeq>
inline void Model<TSeq>::draw(
    DiagramType diagram_type,
    const std::string & fn_output,
    bool self
) {

    ModelDiagram diagram;

    diagram.draw_from_data(
        diagram_type,
        this->get_states(),
        this->get_db().get_transition_probability(false),
        fn_output,
        self
    );

    return;

}

#undef VECT_MATCH
#undef DURCAST
#undef CASES_PAR
#undef CASE_PAR
#undef CHECK_INIT
#endif
