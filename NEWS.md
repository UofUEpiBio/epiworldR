# epiworldR 0.3-2

* Starting version 0.3-0, `epiworldR` is versioned using the same version as the C++ library, `epiworld`.

* Adds the new mixing models `ModelSIRMixing` and `ModelSEIRMixing`.

* Ports the `Entity` class. Entities are used to group agents within a model.

* Refactors `add_tool`, `add_virus`, and `add_entity` simplifying syntax. Now,
  these functions only receive the model and object. Prevalence is
  specified in the object itself. `add_tool_n` and `add_virus_n` are now
  deprecated.

* `globalaction_*` are now defunct. Use `globalevent_*` instead.

* New functions to specify how viruses, tools, and entities are distributed
  among agents: `distribute_viruses`, `distribute_tools`, and `distribute_entities`.


# epiworldR 0.1-0

* Force model to update agents' states when running a simulation.
  This was causing issues when calling `run_multiple()` after a single
  call of `run()`. Reported on [14](https://github.com/UofUEpiBio/epiworldR/issues/14).


# epiworldR 0.0-4

* Added missing checks of tool class when adding a model with `add_too_n`.

* Various small improvements.


# epiworldR 0.0-3

* Added the following models: `ModelSEIRD`, `ModelSEIRDCONN`, `ModelSIRD`, `ModelSIRDCONN`, and `ModelSISD`.

* Fixed a bug reported on issue [6](https://github.com/UofUEpiBio/epiworldR/issues/6).


# epiworldR 0.0-2

* Added a `NEWS.md` file to track changes to the package.

* Fixed bug reported by CRAN about reference to nullptr.

* Renamed arguments across Models in favor of consistency.

* Figures now show virus/tool name instead of id.

* Fixed bug in `run_multiple` and added more tests (C++).

* Redid autoconf and Makevars using RcppArmadillo as a template for checking for OpenMP.
