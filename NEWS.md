# epiworldR 0.2-0 (dev)

* Adds the new mixing models `ModelSIRMixing` and `ModelSEIRMixing`.

* Ports the `Entity` class. Entities are used to group agents within a model.

# epiworldR 0.1-0`

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
