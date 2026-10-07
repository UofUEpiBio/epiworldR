#ifndef EPIWORLD_POSTSAMPLING_MEAT_HPP
#define EPIWORLD_POSTSAMPLING_MEAT_HPP

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::set_post_sampling(PostSamplingFun<TSeq> fun)
{

    post_sampling_fun = std::move(fun);
    post_sampling_on = static_cast< bool >(post_sampling_fun);

    // A model that already has agents sizes the scratch now
    if (post_sampling_on)
        post_sampling_prepare_scratch();

    return *this;

}

/**
 * Sizes the scratch for the current population. Copies, moves and assignments
 * of a model leave the scratch empty, so this is also what makes them safe to
 * continue stepping without a reset.
 */
template<typename TSeq>
inline void Model<TSeq>::post_sampling_prepare_scratch()
{

    if (post_sampling_scratch.counts.size() == population.size())
        return;

    // In 64 bits: on wasm32, size_t is 32 bits and shifting it by 32 is UB
    if (static_cast< uint64_t >(population.size()) >= (uint64_t(1) << 32))
        throw std::length_error(
            "The post-sampling callback supports populations below 2^32 agents."
        );

    // Only `counts` (all zero between steps): the pairs of this step are kept
    post_sampling_scratch.touched.clear();
    post_sampling_scratch.counts.assign(population.size(), 0u);

}

template<typename TSeq>
inline Model<TSeq> & Model<TSeq>::clear_post_sampling()
{
    return set_post_sampling(nullptr);
}

template<typename TSeq>
inline void Model<TSeq>::register_sampled_contact(
    size_t infectious_id,
    size_t contacted_id
)
{
    post_sampling_scratch.pairs.push_back(static_cast< uint32_t >(infectious_id));
    post_sampling_scratch.pairs.push_back(static_cast< uint32_t >(contacted_id));
}

template<typename TSeq>
inline void Model<TSeq>::register_sampled_contacts(
    const size_t * infectious_ids,
    size_t n,
    size_t contacted_id
)
{
    auto & pairs = post_sampling_scratch.pairs;
    for (size_t k = 0u; k < n; ++k)
    {
        pairs.push_back(static_cast< uint32_t >(infectious_ids[k]));
        pairs.push_back(static_cast< uint32_t >(contacted_id));
    }
}

/**
 * Groups the pairs of the step by infectious agent (only the distinct ids are
 * sorted, never an N-sized array) and runs the callback on each group, in
 * ascending id order. The scratch is left clean even if the callback throws.
 */
template<typename TSeq>
inline void Model<TSeq>::post_sampling_dispatch()
{

    auto & sc = post_sampling_scratch;
    const size_t npairs = sc.pairs.size() / 2u;

    if (npairs == 0u)
        return;

    // The callback gets the model: it can clear or replace the callback. The
    // whole dispatch runs on a snapshot, so that takes effect next step.
    const PostSamplingFun<TSeq> fun = post_sampling_fun;

    post_sampling_prepare_scratch();

    try
    {

        // A push visits the carriers in ascending id order, so its pairs are
        // already grouped by infectious agent: no counting, no sorting.
        bool grouped_already = true;
        for (size_t k = 1u; k < npairs; ++k)
            if (sc.pairs[2u * k] < sc.pairs[2u * (k - 1u)])
            {
                grouped_already = false;
                break;
            }

        const size_t nt = [&]() -> size_t {

            sc.grouped.resize(npairs);

            if (grouped_already)
            {

                for (size_t k = 0u; k < npairs; ++k)
                {

                    const uint32_t i = sc.pairs[2u * k];
                    if ((k == 0u) || (i != sc.pairs[2u * (k - 1u)]))
                    {
                        sc.touched.push_back(i);
                        sc.starts.push_back(k);
                    }

                    sc.grouped[k] = sc.pairs[2u * k + 1u];

                }

                sc.starts.push_back(npairs);
                return sc.touched.size();

            }

            // Counting per infectious agent, remembering the distinct ones
            for (size_t k = 0u; k < npairs; ++k)
            {
                const uint32_t i = sc.pairs[2u * k];
                if (sc.counts[i]++ == 0u)
                    sc.touched.push_back(i);
            }

            std::sort(sc.touched.begin(), sc.touched.end());

            // Turning the counts into start positions
            const size_t nd = sc.touched.size();
            sc.starts.resize(nd + 1u);
            size_t pos = 0u;
            for (size_t t = 0u; t < nd; ++t)
            {
                const uint32_t i = sc.touched[t];
                sc.starts[t] = pos;
                pos += sc.counts[i];
                sc.counts[i] = static_cast< uint32_t >(sc.starts[t]);
            }
            sc.starts[nd] = pos;

            // Scattering the contacted ids
            for (size_t k = 0u; k < npairs; ++k)
                sc.grouped[sc.counts[sc.pairs[2u * k]]++] = sc.pairs[2u * k + 1u];

            return nd;

        }();

        // Callbacks must not retain the view: the buffers are reused.
        for (size_t t = 0u; t < nt; ++t)
        {

            const size_t len = sc.starts[t + 1u] - sc.starts[t];
            SampledContactsView view(sc.grouped.data() + sc.starts[t], len);
            fun(&population[sc.touched[t]], view, this);

        }

    }
    catch (...)
    {
        sc.clear();
        throw;
    }

    sc.clear();

}

/**
 * @brief A `PostSamplingFun` that records each batch in the model's
 * `ContactTracing`.
 *
 * @details For every contact, it calls
 * `add_contact(infectious_id, contacted_id, today)`, which is what the
 * built-in models with tracing used to do inline. The model needs
 * `contact_tracing_on()`; turning tracing on does not install this callback.
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline PostSamplingFun<TSeq> make_contact_tracing_post_sampling()
{

    return [](
        Agent<TSeq> * infectious,
        const SampledContactsView & contacts,
        Model<TSeq> * m
    ) -> void {

        if (!m->is_contact_tracing_on())
            return;

        auto & ct = m->get_contact_tracing();
        const size_t id = infectious->get_id();
        const size_t today = static_cast< size_t >(m->today());
        for (const auto c : contacts)
            ct.add_contact(id, c, today);

    };

}

#endif
