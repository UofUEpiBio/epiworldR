#include "../epiworld.hpp"
using namespace epiworld;
#ifndef EPIWORLD_TEMPLATES_VIRUSES_HPP 
#define EPIWORLD_TEMPLATES_VIRUSES_HPP

template<typename TSeq>
inline VirusFun<TSeq> virus_fun_logit(
    std::string name,
    std::vector< size_t > vars,
    std::vector< epiworld_double > coefs
) {

    VirusFun<TSeq> fun_infect = [coefs,vars](
        Agent<TSeq> * agent,
        Virus<TSeq> & virus,
        Model<TSeq> * model
        ) -> VirusFun<TSeq> {

        size_t K = coefs.size();
        epiworld_double res = 0.0;

        #pragma omp simd reduction(+:res)
        for (size_t i = 0u; i < K; ++i)
            res += agent->operator[i] * coefs.at(i);

        return 1.0/(1.0 + std::exp(-res));

    };

    return fun_infect;

}

#endif