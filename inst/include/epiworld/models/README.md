# Models included in epiworld

epiworld features multiple predefined models. The overall structure of the `models`
namespace is:

- It contains factory functions.

- The functions return `epiworld::Model` type objects.

The following models are included:

1. `models::ModelSIS` SIS 

2. `models::ModelSIR` SIR

3. `models::ModelSEIR` SEIR

4. `models::ModelSEIRQuarantine` SEIR with quarantine, isolation and contact tracing

5. `models::ModelSURV` Surveillance

Contacts come from the model's layers (networks, mixing by group, or both), so
there is one model per disease, whatever the contacts.

