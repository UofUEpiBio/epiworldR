#ifndef EPIWORLD_CONTACTMATRIX_BONES_HPP
#define EPIWORLD_CONTACTMATRIX_BONES_HPP

#include <vector>

/**
 * @file contactmatrix-bones.hpp
 * @brief Base class for handling contact matrices in epidemiological
 * models
 * 
 * This class provides functionality for setting and validating contact
 * matrices, which are essential for modeling population mixing in
 * epidemiological simulations. Models that require contact matrices can
 * inherit from this class.
 */
class ContactMatrix 
{
private: 

    std::vector< double > contact_matrix;
    std::vector< double > contact_matrix_backup; ///< Used for resetting the model
    int n_groups = -1;

public:

    ContactMatrix() = default;

    /**
     * @brief Validates the contact matrix size and values
     * @param expected_size The expected size of the contact matrix
     */
    void validate_contact_matrix(size_t expected_size);
    /**
     * @brief Set the contact matrix for population mixing
     * @param cmat Contact matrix specifying interaction rates between groups
     * @param as_backup Whether to use the matrix as a backup (default: true)
     * If set to true, the new contact matrix will be saved as a backup for
     * resetting the model.
     */
    void set_contact_matrix(std::vector< double > cmat, bool as_backup = true);

    /**
     * @brief Get the current contact matrix
     * @return Vector representing the contact matrix
     */
    const std::vector< double > & get_contact_matrix() const;

    /**
     * @brief Get the current contact matrix
     * @return Vector representing the contact matrix
     */
    std::vector< double > & get_contact_matrix_ref();

    /**
     * @brief Get the contact rate between two groups
     * @param i Index of the first group
     * @param j Index of the second group
     * @return Contact rate between group i and group j
     */
    double get_contact_rate(size_t i, size_t j, bool check = true) const;

    /**
     * @brief Get the size of the contact matrix
     * @return Size of the contact matrix
     */
    size_t get_contact_matrix_size() const;
};

#endif