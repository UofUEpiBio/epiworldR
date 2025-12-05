#ifndef EPIWORLD_MODELS_HPP
#define EPIWORLD_MODELS_HPP

/**
 * @defgroup epi_models Epidemiological Models
 * @brief Collection of predefined epidemiological models for disease simulation
 * 
 * This module provides a comprehensive collection of epidemiological models commonly 
 * used in disease simulation and modeling. The models range from basic compartmental 
 * models to specialized scenarios including connected populations, population mixing, 
 * quarantine measures, and disease-specific implementations.
 */

/**
 * @defgroup basic_compartmental Basic Compartmental Models
 * @ingroup epi_models
 * @brief Fundamental compartmental models (SIS, SIR, SEIR)
 * 
 * These are the foundational compartmental models in epidemiology that divide 
 * populations into distinct states based on disease status.
 */

/**
 * @defgroup death_compartmental Models with Death Compartment
 * @ingroup epi_models
 * @brief Models that include a death state (SIRD, SISD, SEIRD)
 * 
 * Extended compartmental models that explicitly track mortality as a separate 
 * compartment in the disease progression.
 */

/**
 * @defgroup connected_models Connected Population Models
 * @ingroup epi_models
 * @brief Models for connected populations (SIRConnected, SEIRConnected, SIRDConnected, SEIRDConnected)
 * 
 * Models designed for populations with explicit network connections between individuals, 
 * allowing for more realistic contact patterns.
 */

/**
 * @defgroup mixing_models Population Mixing Models
 * @ingroup epi_models
 * @brief Models with population mixing (SEIRMixing, SIRMixing)
 * 
 * Models that incorporate heterogeneous mixing patterns between different population 
 * groups or entities, often using contact matrices.
 */

/**
 * @defgroup disease_specific Disease-Specific Models
 * @ingroup epi_models
 * @brief Models tailored for specific diseases (Measles variants, Cholera)
 * 
 * Specialized implementations designed to capture the unique transmission dynamics 
 * and characteristics of specific infectious diseases.
 */

/**
 * @defgroup special_models Specialized Models
 * @ingroup epi_models
 * @brief Other specialized models (DiffNet, SIRLogit, Surveillance)
 * 
 * Models for specific scenarios including diffusion networks, logistic regression-based 
 * transmission, and surveillance systems.
 */

/**
 * @defgroup model_utilities Model Utilities
 * @ingroup epi_models
 * @brief Utility functions for model initialization and global events
 * 
 * Helper functions and global event handlers that can be used across different 
 * epidemiological models for initialization and event management.
 */

/**
 * @namespace epimodels
 * @brief Namespace containing predefined epidemiological models and related functions
 * @ingroup epi_models
 * 
 * This namespace provides a collection of epidemiological models commonly used in
 * disease simulation and modeling. It includes various compartmental models such as
 * SIS, SIR, SEIR, and their variations, as well as specialized models for specific
 * scenarios like connected populations, quarantine measures, and school environments.
 * 
 * The namespace also contains initialization functions, global event handlers, and
 * surveillance mechanisms that can be used across different epidemiological models.
 * 
 * Available model categories include:
 * 
 *   - Basic compartmental models (SIS, SIR, SEIR)
 *   - Models with death compartments (SIRD, SISD, SEIRD)
 *   - Connected population models (SIRConnected, SEIRConnected, SIRDConnected, SEIRDConnected)
 *   - Logistic regression models (SIRLogit)
 *   - Mixing population models (SEIRMixing, SIRMixing)
 *   - Quarantine models (SEIRMixingQuarantine)
 *   - School-specific models (MeaslesSchool)
 *   - Disease-specific models (MeaslesMixing, MeaslesMixingRiskQuarantine)
 *   - Diffusion network models (DiffNet)
 *   - Surveillance systems
 * 
 */
namespace epimodels {

    #include "init-functions.hpp"

    #include "globalevents.hpp"

    #include "sis.hpp"
    #include "sir.hpp"
    #include "seir.hpp"
    #include "surveillance.hpp"
    #include "sirconnected.hpp"
    #include "seirconnected.hpp"
    #include "sird.hpp"
    #include "sisd.hpp"
    #include "seird.hpp"
    #include "sirdconnected.hpp"
    #include "seirdconnected.hpp"
    #include "sirlogit.hpp"
    #include "diffnet.hpp"
    #include "seirmixing.hpp"
    #include "sirmixing.hpp"
    #include "measlesschool.hpp"
    #include "seirmixingquarantine.hpp"
    #include "measlesmixing.hpp"
    #include "measlesmixingriskquarantine.hpp"

}

#endif