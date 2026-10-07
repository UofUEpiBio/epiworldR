#ifndef EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING
#define EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING

/**
 * @brief Functions for sampling viruses
 * 
 */
namespace sampler {

/**
 * @brief Collects the contact probabilities of the neighbors carrying a virus.
 *
 * @details Fills `m->array_double_tmp` and `m->array_virus_tmp` with, for each
 * neighbor of `p` that carries a virus (in the neighbor order), the probability
 * that neighbor's virus is transmitted to `p`, and the virus. Neighbors in a
 * state flagged in `exclude` (if any) are skipped. This is the scan of a pull.
 *
 * During a run, most neighbors do not carry a virus. The model's compact
 * per-agent flag rules them out reading one byte each instead of loading the
 * neighbor's `Agent`; the rest of the loop, including the order in which the
 * arrays are filled, is the same, so the outcome (and the random stream) do not
 * depend on which path runs.
 *
 * @return The number of viruses collected.
 */
template<typename TSeq, bool Record>
inline size_t collect_neighbor_viruses_impl(
    Agent<TSeq> * p,
    Model<TSeq> * m,
    const std::vector< bool > * exclude
)
{

    size_t nviruses_tmp = 0u;

    auto add = [&](Agent<TSeq> * neighbor) -> void {

        // If the state is in the list, exclude it
        if ((exclude != nullptr) && (*exclude)[neighbor->get_state()])
            return;

        auto & v = neighbor->get_virus();
        if (v == nullptr)
            return;

        #ifdef EPI_DEBUG
        if (nviruses_tmp >= m->array_virus_tmp.size())
            throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
        #endif

        // The post-sampling callback sees every eligible contact
        if constexpr (Record)
            m->register_sampled_contact(neighbor->get_id(), p->get_id());

        /* And it is a function of susceptibility_reduction as well */
        m->array_double_tmp[nviruses_tmp] =
            (1.0 - p->get_susceptibility_reduction(v, *m)) *
            v->get_prob_infecting(m) *
            (1.0 - neighbor->get_transmission_reduction(v, *m))
            ;

        m->array_virus_tmp[nviruses_tmp++] = &(*v);

        #ifdef EPI_DEBUG
        if (
            (m->array_double_tmp[nviruses_tmp - 1] < 0.0) |
            (m->array_double_tmp[nviruses_tmp - 1] > 1.0)
            )
        {
            printf_epiworld(
                "[epi-debug] Agent %i's virus has transmission prob outside of [0, 1]: %.4f!\n",
                static_cast<int>(neighbor->get_id()),
                m->array_double_tmp[nviruses_tmp - 1]
                );
        }
        #endif

    };

    auto neighbors = p->neighbors_view(*m);
    const char * carrier = m->get_agent_carrier();

    if (carrier != nullptr)
    {

        const size_t * ids = neighbors.ids();
        const size_t n = neighbors.size();
        for (size_t k = 0u; k < n; ++k)
            if (carrier[ids[k]] != 0)
                add(neighbors.agent(k));

    }
    else
    {

        for (auto * neighbor: neighbors)
            add(neighbor);

    }

    return nviruses_tmp;

}

/**
 * @brief Collects the neighbors' viruses (see `collect_neighbor_viruses_impl`),
 * reporting the contacts to the model's post-sampling callback if there is
 * one. The path is chosen once per call, not per neighbor.
 */
template<typename TSeq>
inline size_t collect_neighbor_viruses(
    Agent<TSeq> * p,
    Model<TSeq> * m,
    const std::vector< bool > * exclude
)
{
    return m->has_post_sampling() ?
        collect_neighbor_viruses_impl<TSeq, true>(p, m, exclude) :
        collect_neighbor_viruses_impl<TSeq, false>(p, m, exclude);
}

/**
 * @brief Update function for susceptible agents that samples from neighbors.
 *
 * @details This is what `make_update_susceptible()` returns. It is a named type
 * (rather than a lambda) so that a `Model` can recognize it among its state
 * update functions and, when it is cheaper, update these agents by pushing
 * infection odds from the agents carrying a virus instead (see
 * `Model::set_transmission_mode()`). Both give the same distribution.
 *
 * Neighbors in any of the states listed in `exclude` are not sources of
 * infection.
 *
 * @tparam TSeq
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
class UpdateSusceptible {
public:

    /// States whose agents cannot transmit (e.g., latent or isolated).
    std::vector< epiworld_fast_uint > exclude;

    explicit UpdateSusceptible(std::vector< epiworld_fast_uint > exclude_ = {})
        : exclude(std::move(exclude_)) {}

    void operator()(Agent<TSeq> * p, Model<TSeq> * m);

private:

    // One entry per state, built on the first call. Held by value: every copy
    // of the function (e.g., in each model copied by run_multiple) builds its
    // own.
    std::vector< bool > exclude_agent_bool;

};

template<typename TSeq>
inline void UpdateSusceptible<TSeq>::operator()(
    Agent<TSeq> * p,
    Model<TSeq> * m
)
{

    // The first time we call it, we need to initialize the vector
    if ((exclude.size() != 0u) && (exclude_agent_bool.size() == 0u))
    {

        exclude_agent_bool.resize(m->get_states().size(), false);
        for (auto s : exclude)
        {
            if (s >= exclude_agent_bool.size())
                throw std::logic_error(
                    std::string("You are trying to exclude a state that is out of range: ") +
                    std::to_string(s) + std::string(". There are only ") +
                    std::to_string(exclude_agent_bool.size()) +
                    std::string(" states in the model.")
                    );

            exclude_agent_bool[s] = true;

        }

    }

    if (p->get_virus() != nullptr)
        throw std::logic_error(
            std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
            std::string("Agent id ") + std::to_string(p->get_id()) +
            std::string(" has a virus.")
            );

    // This computes the prob of getting any neighbor variant
    const size_t nviruses_tmp = collect_neighbor_viruses(
        p, m, exclude.size() == 0u ? nullptr : &exclude_agent_bool
    );

    // No virus to compute
    if (nviruses_tmp == 0u)
        return;

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return;

    p->set_virus(*m, *m->array_virus_tmp[which]);

    return;

}

/**
 * @brief Make a function to sample from neighbors
 * 
 * This is akin to the function default_update_susceptible, with the difference
 * that it will create a function that supports excluding states from the sampling
 * frame. For example, individuals who have acquired a virus can be excluded if
 * in incubation state.
 * 
 * The result is a `sampler::UpdateSusceptible`, which models recognize: they
 * may update these agents by pushing infection odds from the agents carrying
 * a virus instead (same distribution; see `Model::set_transmission_mode()`).
 *
 * @tparam TSeq 
 * @param exclude unsigned vector of states that need to be excluded from the sampling
 * @return The update function.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline std::function<void(Agent<TSeq>*,Model<TSeq>*)> make_update_susceptible(
    std::vector< epiworld_fast_uint > exclude = {}
    )
{

    return UpdateSusceptible<TSeq>(std::move(exclude));

}

/**
 * @brief Make a function to sample from neighbors
 * 
 * This is akin to the function default_update_susceptible, with the difference
 * that it will create a function that supports excluding states from the sampling
 * frame. For example, individuals who have acquired a virus can be excluded if
 * in incubation state.
 * 
 * @tparam TSeq 
 * @param exclude unsigned vector of states that need to be excluded from the sampling
 * @return Virus<TSeq>* of the selected virus. If none selected (or none
 * available,) returns a nullptr;
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> make_sample_virus_neighbors(
    std::vector< epiworld_fast_uint > exclude = {}
)
{
    if (exclude.size() == 0u)
    {

        std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> res = 
            [](Agent<TSeq> * p, Model<TSeq> * m) -> Virus<TSeq>* {

                if (p->get_virus() != nullptr)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has a virus.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nviruses_tmp = 0u;
                for (auto * neighbor: p->neighbors_view(*m)) 
                {
                    
                    if (neighbor->get_virus() == nullptr)
                        continue;

                    auto & v = neighbor->get_virus();

                    #ifdef EPI_DEBUG
                    if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
                        throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                    #endif
                        
                    /* And it is a function of susceptibility_reduction as well */ 
                    m->array_double_tmp[nviruses_tmp] =
                        (1.0 - p->get_susceptibility_reduction(v, *m)) * 
                        v->get_prob_infecting(m) * 
                        (1.0 - neighbor->get_transmission_reduction(v, *m)) 
                        ; 
                
                    m->array_virus_tmp[nviruses_tmp++] = &(*v);
                    
                }

                // No virus to compute
                if (nviruses_tmp == 0u)
                    return nullptr;

                // Running the roulette
                int which = roulette(nviruses_tmp, m);

                if (which < 0)
                    return nullptr;

                return m->array_virus_tmp[which]; 

            };

        return res;


    } else {

        // Making room for the query
        std::shared_ptr<std::vector<bool>> exclude_agent_bool =
            std::make_shared<std::vector<bool>>(0);

        std::shared_ptr<std::vector<epiworld_fast_uint>> exclude_agent_bool_idx =
            std::make_shared<std::vector<epiworld_fast_uint>>(exclude);


        std::function<Virus<TSeq>*(Agent<TSeq>*,Model<TSeq>*)> res = 
            [exclude_agent_bool,exclude_agent_bool_idx](Agent<TSeq> * p, Model<TSeq> * m) -> Virus<TSeq>* {

                // The first time we call it, we need to initialize the vector
                if (exclude_agent_bool->size() == 0u)
                {

                    exclude_agent_bool->resize(m->get_states().size(), false);
                    for (auto s : *exclude_agent_bool_idx)
                    {
                        if (s >= exclude_agent_bool->size())
                            throw std::logic_error(
                                std::string("You are trying to exclude a state that is out of range: ") +
                                std::to_string(s) + std::string(". There are only ") +
                                std::to_string(exclude_agent_bool->size()) + 
                                std::string(" states in the model.")
                                );

                        exclude_agent_bool->operator[](s) = true;

                    }

                }    
                
                if (p->get_virus() != nullptr)
                    throw std::logic_error(
                        std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense! ") +
                        std::string("Agent id ") + std::to_string(p->get_id()) +
                        std::string(" has a virus.")
                        );

                // This computes the prob of getting any neighbor variant
                size_t nviruses_tmp = 0u;
                for (auto * neighbor: p->neighbors_view(*m)) 
                {

                    // If the state is in the list, exclude it
                    if (exclude_agent_bool->operator[](neighbor->get_state()))
                        continue;

                    if (neighbor->get_virus() == nullptr)
                        continue;

                    auto & v = neighbor->get_virus();
                            
                    #ifdef EPI_DEBUG
                    if (nviruses_tmp >= static_cast<int>(m->array_virus_tmp.size()))
                        throw std::logic_error("Trying to add an extra element to a temporal array outside of the range.");
                    #endif
                        
                    /* And it is a function of susceptibility_reduction as well */ 
                    m->array_double_tmp[nviruses_tmp] =
                        (1.0 - p->get_susceptibility_reduction(v, *m)) * 
                        v->get_prob_infecting(m) * 
                        (1.0 - neighbor->get_transmission_reduction(v, *m)) 
                        ; 
                
                    m->array_virus_tmp[nviruses_tmp++] = &(*v);
                    
                }

                // No virus to compute
                if (nviruses_tmp == 0u)
                    return nullptr;

                // Running the roulette
                int which = roulette(nviruses_tmp, m);

                if (which < 0)
                    return nullptr;

                return m->array_virus_tmp[which]; 

            };

        return res;

    }

}

/**
 * @brief Sample from neighbors pool of viruses (at most one)
 * 
 * This function samples at most one virus from the pool of
 * viruses from its neighbors. If no virus is selected, the function
 * returns a `nullptr`, otherwise it returns a pointer to the
 * selected virus.
 * 
 * This can be used to build a new update function (EPI_NEW_UPDATEFUN.)
 * 
 * @tparam TSeq 
 * @param p Pointer to person 
 * @param m Pointer to the model
 * @return Virus<TSeq>* of the selected virus. If none selected (or none
 * available,) returns a nullptr;
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline Virus<TSeq> * sample_virus_single(Agent<TSeq> * p, Model<TSeq> * m)
{

    if (p->get_virus() != nullptr)
        throw std::logic_error(
            std::string("Using the -default_update_susceptible- on agents WITH viruses makes no sense!") +
            std::string("Agent id ") + std::to_string(p->get_id()) +
            std::string(" has a virus.")
            );

    // This computes the prob of getting any neighbor variant
    const size_t nviruses_tmp = collect_neighbor_viruses<TSeq>(p, m, nullptr);

    // No virus to compute
    if (nviruses_tmp == 0u)
        return nullptr;

    #ifdef EPI_DEBUG
    m->get_db().n_transmissions_potential++;
    #endif

    // Running the roulette
    int which = roulette(nviruses_tmp, m);

    if (which < 0)
        return nullptr;

    #ifdef EPI_DEBUG
    m->get_db().n_transmissions_today++;
    #endif

    return m->array_virus_tmp[which]; 
    
}

}

#endif