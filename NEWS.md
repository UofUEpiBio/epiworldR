# epiworldR 0.11.3.0

* Updates to `epiworld` version 0.11.3, including a patch to the `ModelSEIRMixingQuarantine` model. It was underestimating the outbreak size. This is a similar issue detected in the collection of the measles models.


# epiworldR 0.11.2.0

* Updates to `epiworld` version 0.11.2 including minor patches to avoid warnings and notes from clang.

## Bug Fixes

These changes will be reflected in the `measles` R package (so no changes for the `epiworldR` package):

* The measles school model was not including the number of recovered agents in the population that it could sample from, which resulted in a minor overestimation of the outbreak size (about 5 to 10 percent).

* The measles mixing models were assuming the number of available agents did not change when sampling contacts. This has been corrected. Outbreak sizes are now expected to be larger on average.


# epiworldR 0.11.0.1

* Minor release fixing errors in tests in some OSs (reported by CRAN).


# epiworldR 0.11.0.0

## New Features

* New distribution functions `distribute_virus_to_entities()` and `distribute_tool_to_entities()` that allow distributing viruses and tools to agents based on their entity membership (new in epiworld 0.11.0).

* New function `add_entities_from_dataframe()` that allows adding multiple entities to a model based on a data.frame (wrapper around `add_entity()`).

* New functions `get_outbreak_size()` and `get_active_cases()` allow extracting timeseries of the outbreak size and number of active cases (individuals with a virus) (new in epiworld 0.11.0).

* New function `get_hospitalizations()` that allows extracting timeseries of hospitalizations from the model results (new in epiworld 0.11.0).

## Bug Fixes

* The saver was not correctly saving the history of tools and viruses when using `virus_hist` and `tools_hist`, respectively. This has been fixed (new in epiworld 0.11.0).

## Breaking Changes

* **Measles models have been removed from epiworldR**. The measles-specific models (`ModelMeaslesSchool`, `ModelMeaslesMixing`, `ModelMeaslesMixingRiskQuarantine`, and the deprecated `ModelMeaslesQuarantine`) have been moved to a separate package. Users who need measles models should install the [`measles`](https://github.com/UofUEpiBio/measles) package.


# epiworldR 0.10.0.0

## New Features

* New model `ModelMeaslesMixingRiskQuarantine` that implements a measles model with mixing, risk-based contacts, and quarantine has been added. This is part of a joint effort with our partners at the Utah Department of Health and Human Services (DHHS).

* Added `print` and `plot` methods for `epiworld_diagram` objects. This allows easy visualization of model diagrams generated with `DiagrammeR`.

## Improvements

* The function `run_multiple_get_results()` has two new arguments: `freader` and `...`, which allows the user to specify a custom function to read the results of each simulation (e.g., `data.table::fread`).

## Bug Fixes

* The rewiring algorithm was not changing the edges, but only the weights of the edges. This has been fixed in epiworld 0.9.4.

* The saver associated with `tools_hist` had a bug that caused it to not save the history of tools correctly. This has been fixed in epiworld 0.10.0. This was detected during work with @margowheatley and @AliseMendoza.

## Misc

* The `ModelMeaslesQuarantine` model has been updated to warn the user if the deprecated parameter `vax_improved_recovery` is used. This parameter is no longer used and has been removed.

* Various legacy C++ files have been removed from the `inst/include/epiworld` folder.

* We now include model diagrams in the documentation of each model (new feature of `epiworld` 0.10.0).

* Various improvements to the documentation of the package and models.


# epiworldR 0.9.1.0

## Improvements

* The C++ library has been updated to version 0.9.1, which includes several speed improvements. Some preliminary tests show x10 speedups.

* New model `ModelSEIRMixingQuarantine` that implements a SEIR model with quarantine and mixing has been added.

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
