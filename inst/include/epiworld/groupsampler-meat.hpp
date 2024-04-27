#ifndef GROUPSAMPLER_MEAT_HPP
#define GROUPSAMPLER_MEAT_HPP

template<typename TSeq>
inline GroupSampler<TSeq>::GroupSampler(
        const std::vector< double > & contact_matrix_,
        const std::vector< size_t > & group_sizes_,
        bool normalize
    ): contact_matrix(contact_matrix_), group_sizes(group_sizes_) {


        this->cumulate.resize(contact_matrix.size());
        std::fill(cumulate.begin(), cumulate.end(), 0.0);

        // Cumulative sum
        for (size_t j = 0; j < group_sizes.size(); ++j)
        {
            for (size_t i = 0; i < group_sizes.size(); ++i)
                cumulate[idx(i, j, true)] += 
                    cumulate[idx(i, j - 1, true)] +
                    contact_matrix[idx(i, j)];
        }

        if (normalize)
        {
            for (size_t i = 0; i < group_sizes.size(); ++i)
            {
                double sum = 0.0;
                for (size_t j = 0; j < group_sizes.size(); ++j)
                    sum += contact_matrix[idx(i, j, true)];
                for (size_t j = 0; j < group_sizes.size(); ++j)
                    contact_matrix[idx(i, j, true)] /= sum;
            }
        }

    };

template<typename TSeq>
int GroupSampler<TSeq>::sample_1(
    Model<TSeq> * model,
    const int origin_group
    )
{

    // Random number
    double r = model->runif();

    // Finding the group
    size_t j = 0;
    while (r > cumulate[idx(origin_group, j, true)])
        ++j;

    // Adjusting the prob
    r = r - (j == 0 ? 0.0 : cumulate[idx(origin_group, j - 1, true)]);

    int res = static_cast<int>(
        std::floor(r * group_sizes[j])
    );

    // Making sure we are not picking outside of the group
    if (res >= static_cast<int>(group_sizes[j]))
        res = static_cast<int>(group_sizes[j]) - 1;

    return model->get_entities()[j][res]->get_id();

}

template<typename TSeq>
void GroupSampler<TSeq>::sample_n(
    Model<TSeq> * model,
    std::vector< int > & sample,
    const int origin_group,
    const int nsamples
)
{

    for (int i = 0; i < nsamples; ++i)
        sample[i] = sample_1(model, origin_group);

    return;

}

#endif