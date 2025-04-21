# epiworldR 0.8.1.0

* Added the `ModelMeaslesQuarantine`.

* Fixes an error in the index case calculation of the Rt function. Getting the Rt was not recovering the index cases with no transmissions.

* Adds the capability of drawing a model diagram using `mermaid`.

* Added an advanced usage vignette.

* Adding input parameter checks (more extensive).


# epiworldR 0.6.1.0

* Updates to reflect changes in the `epiworld` C++ library (mostly bug fixes)

* Package now requires R version >=4.1.0, because it uses the pipe `|>`


# epiworldR 0.6.0.0

* The package now includes the `LFMCMC` module that implements
  the likelihood-free Markov Chain Monte Carlo algorithm. This
  module is used to estimate the parameters of the models.

* The new function `add_param()` allows the user to add parameters
  to the model.

* The new function `rm_globalevent()` allows the user to remove
  global events from the model.

* The function `today()` returns the current day (step) of the
  simulation.

* We changed the versioning system. To allow the R package to increase
  version number while preserving epiworld (C++) versioning, we added a fourth
  number that indicates R-only patches (similar to RcppArmadillo).


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
