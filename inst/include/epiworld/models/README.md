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

A model's states describe the course of infection; its layers (networks,
mixing by group, or both) describe who meets whom. So there is one model per
compartmental structure, whatever the contacts, and a model can carry several
viruses (`add_virus()`), which share its states.

