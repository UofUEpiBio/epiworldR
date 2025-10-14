# epiworldR 0.9.4.0

## Bug Fixes

* The rewiring algorithm was not changing the edges, but only the weights of the edges. This has been fixed in epiworld 0.9.4.


# epiworldR 0.9.1.0

## Improvements

* The C++ library has been updated to version 0.9.1, which includes several speed improvements. Some preliminary tests show x10 speedups.

* New model `ModelSEIRMixingQuarantine` that implements a SEIR model with quarantine and mixing has been added.

* New model `ModelMeasles` that implements a measles based on the `ModelSEIRMixingQuarantine` model has been added.

## Other changes

* The Model `ModelMeaslesQuarantine` has been deprecated in favor of `ModelMeaslesSchool` (change of name).


# epiworldR 0.8.3.0 (CRAN)

## Improvements

* epiworld was storing empty entries in the transition matrix. This was making reading the saved data significantly slower.

## Bug Fixes

* The `run_multiple()` function was failing to correctly backup the agents' edgelist in the case of network models. This was due to a bug in the C++ library that was not correctly copying the edgelist. This has been fixed.

* When multiple transitions happened in a single step, e.g., I->E->S, the model was only recording E->S, but not I->S (which is the correct).


# epiworldR 0.8.2.0

## Changes

* In the `ModelMeaslesQuarantine` model, isolated individuals now have a different parameter for the number of days in isolation. Previously, the same parameter was used for both isolated and quarantined individuals. The new parameter is `isolation_period`.

* The `ModelMeaslesQuarantine` model changed a couple of parameter names to improve consistency. We also improved the documentation of the model.

* The mixing models were using the wrong sampling scheme for the contact event between agents. The new implementation allows using the regular basic reproductive number calculation for its calibration.

* The functions `distribute_(virus|tool)_to_set` now have an additional argument that allows restricting the set of agents.

## Other Changes

* Added more tests to the `ModelMeaslesQuarantine` model, as well as to the `set_distribution*`, `distribute_*`, and `get_agents_states` functions.

* A large overhaul memory management was performed (C++ side), decreasing the memory footprint by one order of magnitude in some cases.


# epiworldR 0.8.1.0

## New Features

* Adds the `ModelMeaslesQuarantine` model, which implements a
  measles model with quarantine.

* Adds the capability of generating model diagrams using `mermaid` through the `draw_mermaid()` model function and the `ModelDiagram` set of functions.

## Enhancements

* Adds an advanced usage vignette.

* Adds more extensive input parameter checks during model creation.

## Bug Fixes

* Fixes an error in the index case calculation of the Rt function. Getting the Rt was not recovering the index cases with no transmissions.

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
