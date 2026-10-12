args <- commandArgs(TRUE)
library(epiworldR, lib.loc = args[1])
source("inst/tinytest/helpers/network-recipes.R")
variant <- args[2]
round <- as.integer(args[3])
out <- args[4]
kinds <- c("SIR", "SIS", "SEIR", "SIRD", "SISD", "SEIRD")
timings <- list()
snapshots <- list()
record <- function(kind, n, mode, phase, seconds) {
  timings[[length(timings) + 1L]] <<- data.frame(
    kind, n, mode, phase, variant, round, seconds)
}
for (kind in kinds) {
  constructor <- get(paste0("Model", kind))
  params <- recipe_args(kind)
  for (j in seq_len(20L)) invisible(do.call(constructor, params))
  gc()
  elapsed <- system.time(for (j in seq_len(1000L))
    invisible(do.call(constructor, params)))[["elapsed"]]
  record(kind, 0L, "none", "construction", elapsed / 1000)
  for (n in c(1000L, 10000L)) for (mode in c("pull", "auto")) {
    m <- recipe_model(kind, n = n, mode = mode)
    run(m, ndays = 100L, seed = 1L) # warm up; graph creation is outside timing
    gc()
    repetitions <- if (n == 1000L) 150L else 25L
    elapsed <- system.time(for (j in seq_len(repetitions))
      run(m, ndays = 100L, seed = 1000L * round + j))[["elapsed"]]
    record(kind, n, mode, "simulation", elapsed / repetitions)
    # Final repetition: exact total, virus, transition and transmission data.
    snapshots[[paste(kind, n, mode)]] <- recipe_snapshot(m)
  }
}
write.csv(do.call(rbind, timings),
  file.path(out, paste0("timings-", variant, "-", round, ".csv")), row.names = FALSE)
saveRDS(snapshots, file.path(out, paste0("snapshots-", variant, "-", round, ".rds")))
if (round == 1L) capture.output(sessionInfo(),
  file = file.path(out, paste0("session-", variant, ".txt")))
cat(variant, "round", round, "complete\n")
