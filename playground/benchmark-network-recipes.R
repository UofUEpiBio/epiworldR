# Run inside the devcontainer after clean installs of base and PR into two
# separate libraries. See benchmark-network-recipes.md for exact commands.
args <- commandArgs(TRUE)
if (length(args) != 3L)
  stop("Usage: Rscript playground/benchmark-network-recipes.R BASE_LIB PR_LIB OUTPUT_DIR")
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[3])
worker <- normalizePath("playground/benchmark-network-recipes-worker.R")
for (round in seq_len(9L)) {
  # Alternate which build goes first; both use the same seeds in each round.
  for (variant in if (round %% 2L) c("cpp", "recipe") else c("recipe", "cpp")) {
    lib <- normalizePath(args[if (variant == "cpp") 1L else 2L])
    status <- system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(worker, lib, variant, as.character(round), out)))
    if (status != 0L) stop("Benchmark worker failed")
  }
}
files <- list.files(out, pattern = "^timings-.*csv$", full.names = TRUE)
timings <- do.call(rbind, lapply(files, read.csv))
write.csv(timings, file.path(out, "timings.csv"), row.names = FALSE)
summary <- aggregate(seconds ~ kind + n + mode + phase + variant, timings, median)
summary <- reshape(summary, direction = "wide", timevar = "variant",
  idvar = c("kind", "n", "mode", "phase"))
summary$ratio <- summary$seconds.recipe / summary$seconds.cpp
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
for (round in seq_len(9L)) {
  a <- readRDS(file.path(out, paste0("snapshots-cpp-", round, ".rds")))
  b <- readRDS(file.path(out, paste0("snapshots-recipe-", round, ".rds")))
  stopifnot(identical(a, b))
}
cat("All paired benchmark trajectories are identical.\n")
