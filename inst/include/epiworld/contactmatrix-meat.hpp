#ifndef EPIWORLD_CONTACTMATRIX_MEAT_HPP
#define EPIWORLD_CONTACTMATRIX_MEAT_HPP

#include "contactmatrix-bones.hpp"
#include <stdexcept>

inline void ContactMatrix::validate_contact_matrix(size_t expected_size)
{
    if (contact_matrix.size() != expected_size * expected_size)
        throw std::length_error(
            std::string("Contact matrix size is ") +
            std::to_string(contact_matrix.size()) +
            std::string(", but expected size is ") +
            std::to_string(expected_size * expected_size) + "."
        );

    for (int i = 0; i < n_groups; ++i)
    {
        for (int j = 0; j < n_groups; ++j)
        {
            if (get_contact_rate(i, j, false) < 0.0)
                throw std::range_error(
                    std::string("The contact matrix must be non-negative. ") +
                    std::to_string(this->get_contact_rate(i, j, false)) +
                    std::string(" < 0.")
                    );
        }
    }
}

inline void ContactMatrix::set_contact_matrix(
    std::vector< double > cmat
)
{
    n_groups = static_cast<int>(std::sqrt(cmat.size()));
    if (n_groups * n_groups != static_cast<int>(cmat.size()))
        throw std::invalid_argument(
            "Contact matrix size is not a perfect square, cannot determine number of groups."
        );

    contact_matrix = cmat;
    return;
};

inline const std::vector< double > & ContactMatrix::get_contact_matrix() const
{
    return contact_matrix;
}

inline std::vector< double > & ContactMatrix::get_contact_matrix_ref()
{
    return contact_matrix;
}

inline double ContactMatrix::get_contact_rate(
    size_t i, size_t j, bool check
) const
{
    if (check && (
        (static_cast<int>(i) >= n_groups) ||
        (static_cast<int>(j) >= n_groups)
    ))
        throw std::out_of_range(
            std::string("Group indices out of range. ") +
            std::to_string(i) + ", " + std::to_string(j) +
            std::string(" >= ") + std::to_string(n_groups)
        );
    return contact_matrix[j * n_groups + i];
}

inline size_t ContactMatrix::get_contact_matrix_size() const
{
    return contact_matrix.size();
}

#endif