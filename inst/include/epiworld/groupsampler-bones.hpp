#ifndef GROUPSAMPLER_BONES_HPP
#define GROUPSAMPLER_BONES_HPP

/**
 * @brief Weighted sampling of groups
 */
template<typename TSeq> 
class GroupSampler {

private:
    
    std::vector< double > contact_matrix; ///< Contact matrix between groups
    std::vector< size_t > group_sizes;    ///< Sizes of the groups
    std::vector< double > cumulate;       ///< Cumulative sum of the contact matrix (row-major for faster access)

    /**
     * @brief Get the index of the contact matrix
     * 
     * The matrix is a vector stored in column-major order.
     * 
     * @param i Index of the row
     * @param j Index of the column
     * @return Index of the contact matrix
     */
    inline int idx(const int i, const int j, bool rowmajor = false) const
    {
        
        if (rowmajor)
            return i * group_sizes.size() + j;
        
        return j * group_sizes.size() + i; 

    }

public:

    GroupSampler() {};

    GroupSampler(
        const std::vector< double > & contact_matrix_,
        const std::vector< size_t > & group_sizes_,
        bool normalize = true
    );

    int sample_1(
        Model<TSeq> * model,
        const int origin_group
        );

    void sample_n(
        Model<TSeq> * model,
        std::vector< int > & sample,
        const int origin_group,
        const int nsamples
    );

};

#endif