# Network model recipes versus C++ constructors

The baseline is epiworldR commit `ce1cb0f0` (the original compiled constructors).
The recipe build replaces `ModelSIR`, `ModelSIS`, `ModelSEIR`, `ModelSIRD`,
`ModelSISD`, and `ModelSEIRD` with R assembly over native building blocks.
Connected, mixing, surveillance, diffusion, and logit models are outside this
change. All six public constructors retain their signatures and S3 classes.

## Reproduce

Use this repository's `.devcontainer/Containerfile` and `devcontainer.json`:

```sh
devcontainer up --workspace-folder . --docker-path podman
devcontainer exec --workspace-folder . --docker-path podman bash
```

Inside the devcontainer, from the PR checkout (normally `/workspaces/epiworldR`):

```sh
mkdir -p /tmp/epiworldR-baseline /tmp/epiworldR-lib-baseline /tmp/epiworldR-lib-recipes
git archive ce1cb0f0 | tar -x -C /tmp/epiworldR-baseline
(cd /tmp/epiworldR-baseline && MAKEFLAGS=-j4 R CMD INSTALL --preclean --clean --no-docs --library=/tmp/epiworldR-lib-baseline .)
MAKEFLAGS=-j4 R CMD INSTALL --preclean --clean --no-docs --library=/tmp/epiworldR-lib-recipes .
Rscript playground/benchmark-network-recipes.R /tmp/epiworldR-lib-baseline /tmp/epiworldR-lib-recipes /tmp/network-recipe-benchmark
```

The driver alternates the baseline/recipe order across nine rounds. Each
worker is a fresh R process loading one installed package. Construction times
are per constructor call over batches of 1,000, after 20 warmup calls.
Simulation times are per 100-day run over 150 seeded repetitions (1,000 agents)
or 25 repetitions (10,000 agents), after a
warmup run, on a fixed degree-10 undirected ring of 1,000 or 10,000 agents.
Graphs and models are created outside the simulation timer. Both `pull` and
`auto` transmission modes are measured; each round uses the same seeds in
both builds. Timing uses `system.time()[["elapsed"]]`, includes run/reset
overhead, and excludes R process startup and loading. Reported times are
medians across nine rounds; small differences should not be generalized
beyond these workloads or this container/VM.

The driver also checks exact equality of total, virus, transition, and
transmission histories from the final timed repetition of each configuration
(216 paired trajectories). Separate regression fingerprints cover 72 combinations
of six models, three transmission modes, queue on/off, and default/custom
initialization. Custom initialization includes parameter edits and treatments.
Regenerate that fixture using the baseline library only:

```sh
R_LIBS=/tmp/epiworldR-lib-baseline Rscript -e 'library(epiworldR); source("inst/tinytest/helpers/network-recipes.R"); reference <- recipe_cases; reference$fingerprint <- vapply(recipe_reference(), recipe_fingerprint, character(1)); write.csv(reference, "inst/tinytest/helpers/network-recipes-reference.csv", row.names = FALSE)'
```

## Compatibility details

`update_fun_infected()` uses the core recovery/death update, including tool
effects and virus removal. `update_fun_incubation()` keeps the virus during
the state transition. `update_fun_recovery()` preserves SEIR's direct
parameter-based recovery, which does not apply treatment recovery effects.
Susceptible updates with no exclusions use the core default function pointer
so the transmission dispatcher recognizes the same optimized path.

One configurable native builder type retains the existing initial-state
factories and supports cloning and parallel simulation. No six model-specific
C++ constructors remain in the package's registered bindings.

This refactor deliberately preserves the original SISD model's hard-coded
0.01 death probability (its `Death rate` parameter currently has no effect)
and the existing SEIR/SEIRD transmission from exposed carriers. Changing
either would require a separate behavior change and invalidate exact
trajectory comparisons.

## Results

Measured in the repository's devcontainer on Linux aarch64 in the macOS Podman
VM (10 virtual CPUs), R 4.6.1, GCC 13.3.0, default release `-O2` flags with
OpenMP enabled. The simulation benchmark uses single-threaded `run()`. No
compilation or test suite ran concurrently with the timed measurements.
Raw timings, all 24 simulation configurations, and session information are in
[`benchmark-network-recipes-results/`](benchmark-network-recipes-results/).

For 10,000 agents, 100 days, default `auto` transmission (median milliseconds):

| Model | C++ constructor | R recipe | Recipe / C++ |
| --- | ---: | ---: | ---: |
| SIR | 6.68 | 6.76 | 1.012 |
| SIS | 92.40 | 93.48 | 1.012 |
| SEIR | 7.84 | 7.72 | 0.985 |
| SIRD | 7.04 | 7.28 | 1.034 |
| SISD | 74.20 | 75.56 | 1.018 |
| SEIRD | 9.16 | 8.84 | 0.965 |

Across both sizes and both modes, ratios range from 0.965 to 1.078. The largest
observed slowdown is SEIRD with 1,000 agents in `pull` mode: 0.853 versus
0.920 ms (+7.8%). These measurements support similar simulation performance,
not a claim of zero overhead or statistically established equivalence.

Construction is more expensive because the recipe makes multiple R/native
calls (median microseconds per constructor, network creation excluded):

| Model | C++ constructor | R recipe | Extra time |
| --- | ---: | ---: | ---: |
| SIR | 8 | 55 | 47 |
| SIS | 8 | 52 | 44 |
| SEIR | 9 | 76 | 67 |
| SIRD | 9 | 59 | 50 |
| SISD | 8 | 57 | 49 |
| SEIRD | 12 | 81 | 69 |

The stripped shared library shrinks from **1,318,872 to 1,253,240 bytes**:
**65,632 bytes (64.1 KiB, 4.98%) smaller**. Unstripped sizes are 25,460,336 and
25,366,040 bytes (94,296 bytes smaller). These are shared-library sizes, not
installed-package or archive sizes. Both builds used the same toolchain and
flags. To reproduce the stripped comparison without changing the installed
libraries:

```sh
strip --strip-unneeded -o /tmp/baseline-stripped.so /tmp/epiworldR-lib-baseline/epiworldR/libs/epiworldR.so
strip --strip-unneeded -o /tmp/recipes-stripped.so /tmp/epiworldR-lib-recipes/epiworldR/libs/epiworldR.so
stat -c '%s %n' /tmp/baseline-stripped.so /tmp/recipes-stripped.so
```

All 216 paired benchmark snapshots match exactly. The 72 regression cases
also match the original C++ constructors, including the preserved SISD death
probability behavior. The fixtures store compact fingerprints of full
version-2 serialized snapshots to avoid shipping large simulation datasets.

Validation in the same devcontainer: all **888 tinytest results passed**,
including cloning, changing initial states after cloning, and parallel runs.
`R CMD check --no-manual --no-build-vignettes --no-tests` completed with no
errors and two warnings about omitted vignette HTML/PDF outputs; all examples
and extracted vignette R code passed. Tests were run separately before check.
