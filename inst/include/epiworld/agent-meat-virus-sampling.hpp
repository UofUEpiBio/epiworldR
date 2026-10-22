#ifndef EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING
#define EPIWORLD_AGENT_MEAT_VIRUS_SAMPLING

/**
 * @brief Functions for sampling viruses
 * 
 */
namespace sampler {

/**
 * @brief Update function for susceptible agents that samples from their
 * neighbors in every network layer (see `Model::sample_virus()`).
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

    Virus<TSeq> * virus = m->sample_virus(
        p, exclude.size() == 0u ? nullptr : &exclude_agent_bool
    );

    if (virus != nullptr)
        p->set_virus(*m, *virus);

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

}

#endif