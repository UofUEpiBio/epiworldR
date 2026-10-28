#ifndef EPIWORLD_MODELS_HPP
#define EPIWORLD_MODELS_HPP

/**
 * @defgroup epi_models Epidemiological Models
 * @brief Collection of predefined epidemiological models for disease simulation
 * 
 * This module provides a comprehensive collection of epidemiological models commonly
 * used in disease simulation and modeling. A model is a disease (its states,
 * parameters and update rules); where contacts come from is up to its layers
 * (networks, mixing by group, or both; see `Model::add_layer()`). The models
 * range from basic compartmental models to quarantine measures and
 * disease-specific implementations.
 * 
 * **Implementation details:**
 * <a href="../impl/core-simulation-models.md">Core Simulation Models</a>,
 * <a href="../impl/extending-the-library-with-models.md">Extending the Library with Models</a>
 */

/**
 * @defgroup basic_compartmental Basic Compartmental Models
 * @ingroup epi_models
 * @brief Fundamental compartmental models (SIS, SIR, SEIR, SEIR with quarantine)
 * 
 * These are the foundational compartmental models in epidemiology that divide 
 * populations into distinct states based on disease status.
 * 
 * **Implementation details:**
 * <a href="../impl/core-simulation-models.md">Core Simulation Models</a>
 */

/**
 * @defgroup death_compartmental Models with Death Compartment
 * @ingroup epi_models
 * @brief Models that include a death state (SIRD, SISD, SEIRD)
 * 
 * Extended compartmental models that explicitly track mortality as a separate 
 * compartment in the disease progression.
 * 
 * **Implementation details:**
 * <a href="../impl/core-simulation-models.md">Core Simulation Models</a>
 */

/**
 * @defgroup disease_specific Disease-Specific Models
 * @ingroup epi_models
 * @brief Models tailored for specific diseases (Measles variants, Cholera)
 * 
 * Specialized implementations designed to capture the unique transmission dynamics 
 * and characteristics of specific infectious diseases.
 * 
 * **Implementation details:**
 * <a href="../impl/core-simulation-models.md">Core Simulation Models</a>,
 * <a href="../impl/quarantine-isolation-and-contact-tracing.md">Quarantine, Isolation, and Contact Tracing</a>
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
 * 
 * **Implementation details:**
 * <a href="../impl/events-and-multi-event-handling.md">Events and Multi-Event Handling</a>,
 * <a href="../impl/extending-the-library.md">Extending the Library</a>
 */

/**
 * @namespace epimodels
 * @brief Namespace containing predefined epidemiological models and related functions
 * @ingroup epi_models
 * 
 * This namespace provides a collection of epidemiological models commonly used in
 * disease simulation and modeling. It includes various compartmental models such as
 * SIS, SIR, SEIR, and their variations, as well as specialized models for specific
 * scenarios like quarantine measures and school environments.
 * 
 * The namespace also contains initialization functions, global event handlers, and
 * surveillance mechanisms that can be used across different epidemiological models.
 * 
 * Available model categories include:
 * 
 *   - Basic compartmental models (SIS, SIR, SEIR)
 *   - Models with death compartments (SIRD, SISD, SEIRD)
 *   - Logistic regression models (SIRLogit)
 *   - Quarantine models (SEIRQuarantine)
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
    #include "sird.hpp"
    #include "sisd.hpp"
    #include "seird.hpp"
    #include "sirlogit.hpp"
    #include "diffnet.hpp"
    #include "seirquarantine.hpp"
}

#endif