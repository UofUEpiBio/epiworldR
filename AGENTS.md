# Instructions for coding agents

`epiworldR` is an R package that uses cpp11 to wrap the C++ library
[epiworld](https://github.com/UofUEpiBio/epiworld), whose headers live in
`inst/include/epiworld`. See `DEVELOPMENT.md` for versioning and release notes.

## Use the devcontainer for everything that runs R or compiles

Compile, install, test, run roxygen2/cpp11 registration, render vignettes and
run `R CMD check` **inside the devcontainer** (`.devcontainer/`), not with the
host's R. The container pins the R, roxygen2, cpp11, tinytest and quarto
versions the project expects; host installs drift (e.g., an older tinytest
without `expect_length()`, or a different roxygen2 that rewrites every `.Rd`).

Start it once, then run commands through `devcontainer exec`. The repository is
mounted at `/workspaces/epiworldR`.

```bash
devcontainer up --workspace-folder .
devcontainer exec --workspace-folder . bash -c 'cd /workspaces/epiworldR && <command>'
```

With Podman instead of Docker, add `--docker-path podman` to both calls.

Common commands (run inside the container):

```bash
# Regenerate cpp11 bindings (after changing a [[cpp11::register]] function)
# and the documentation
Rscript -e 'cpp11::cpp_register(); roxygen2::roxygenise()'

# Clean install
rm -f src/*.o src/*.so && R CMD INSTALL --preclean --no-docs .

# Tests
Rscript -e 'tinytest::test_package("epiworldR", testdir = "inst/tinytest")'

# One vignette
quarto render vignettes/bubbles.qmd --to html
```

Always install from a clean `src/`. Header-only changes (e.g., syncing
`inst/include/epiworld`) are not tracked by R's build, so stale object files get
linked in silently. Since the repository is bind-mounted, objects compiled on
the host and in the container also end up in the same `src/`. Remove
`src/*.o src/*.so` when you are done.

## Syncing the C++ library

`make update` copies `include/epiworld` from the upstream `master` branch (use
`EPIWORLD_BRANCH=<tag>` for a release tag). The first three version numbers in
`DESCRIPTION` follow the C++ library, and the fourth is for R-only changes (e.g.,
epiworld 0.17.0 -> `0.17.0-0`). Summarize user-visible changes and bug fixes in
`NEWS.md`, describing them in terms of the R API.
