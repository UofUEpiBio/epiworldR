#ifndef EPIWORLD_SEQ_PROCESSING_HPP 
#define EPIWORLD_SEQ_PROCESSING_HPP

/**
 * @brief Hasher function to turn the sequence into an integer vector
 * 
 * @tparam TSeq 
 * @param x 
 * @return std::vector<int> 
 */
template<typename TSeq>
inline std::vector<int> default_seq_hasher(const TSeq & x);

template<>
inline std::vector<int> default_seq_hasher<std::vector<int>>(const std::vector<int> & x) {
    return x;
}

template<>
inline std::vector<int> default_seq_hasher<std::vector<bool>>(const std::vector<bool> & x) {
    std::vector<int> ans(x.size());
    size_t j = 0;
    for (const auto & i : x)
        ans[j++] = i? 1 : 0;
    return ans;
}

template<>
inline std::vector<int> default_seq_hasher<int>(const int & x) {
    return {x};
}

template<>
inline std::vector<int> default_seq_hasher<bool>(const bool & x) {
    return {x ? 1 : 0};
}

/**
 * @brief Default way to write sequences
 * 
 * @tparam TSeq 
 * @param seq 
 * @return std::string 
 */
template<typename TSeq = EPI_DEFAULT_TSEQ>
inline std::string default_seq_writer(const TSeq & seq);

template<>
inline std::string default_seq_writer<std::vector<int>>(
    const std::vector<int> & seq
) {

    std::string out = "";
    for (const auto & s : seq)
        out = out + std::to_string(s);

    return out;

}

template<>
inline std::string default_seq_writer<std::vector<bool>>(
    const std::vector<bool> & seq
) {

    std::string out = "";
    for (const auto & s : seq)
        out = out + (s ? "1" : "0");

    return out;

}

template<>
inline std::string default_seq_writer<bool>(
    const bool & seq
) {

    return seq ? "1" : "0";

}

template<>
inline std::string default_seq_writer<int>(
    const int & seq
) {

    return std::to_string(seq);

}



#endif