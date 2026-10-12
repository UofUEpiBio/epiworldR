# Shared by the regression test and benchmark; no timing includes this setup.
recipe_args <- function(kind) {
  args <- list(name = "virus", prevalence = 0.1,
    transmission_rate = 0.12, recovery_rate = 0.2)
  if (kind %in% c("SEIR", "SEIRD")) args$incubation_days <- 3
  if (kind %in% c("SIRD", "SISD", "SEIRD")) args$death_rate <- 0.04
  args
}

recipe_initial <- list(SIR = 0.2, SIRD = c(0.2, 0.1),
  SEIR = c(0.4, 0.2), SEIRD = c(0.4, 0.2, 0.1))

recipe_model <- function(kind, n = 120L, mode = "auto", custom = FALSE,
                         queue = TRUE) {
  m <- do.call(get(paste0("Model", kind)), recipe_args(kind))
  # A deterministic ring with degree 10: no network RNG or cross-version
  # graph generation differences can confound the trajectory comparisons.
  from <- rep(seq_len(n) - 1L, each = 5L)
  to <- (from + rep(seq_len(5L), n)) %% n
  agents_from_edgelist(m, source = from, target = to, size = n, directed = FALSE)
  verbose_off(m)
  set_transmission_mode(m, mode)
  if (!queue) queuing_off(m)
  if (custom) {
    if (!is.null(recipe_initial[[kind]]))
      initial_states(m, recipe_initial[[kind]])
    set_param(m, "Transmission rate", 0.15)
    set_param(m, "Recovery rate", 0.25)
    if (kind %in% c("SEIR", "SEIRD")) set_param(m, "Incubation days", 2)
    if (kind %in% c("SIRD", "SISD", "SEIRD")) set_param(m, "Death rate", 0.08)
    add_tool(m, tool("treatment", prevalence = 0.3, as_proportion = TRUE,
      susceptibility_reduction = 0.2, transmission_reduction = 0.1,
      recovery_enhancer = 0.3, death_reduction = 0.4))
  }
  m
}

recipe_snapshot <- function(m) {
  list(states = get_states(m), total = get_hist_total(m),
    virus = get_hist_virus(m), transitions = get_hist_transition_matrix(m),
    transmissions = get_transmissions(m))
}

recipe_cases <- expand.grid(kind = c("SIR", "SIS", "SEIR", "SIRD", "SISD", "SEIRD"),
  mode = c("pull", "push", "auto"), custom = c(FALSE, TRUE),
  queue = c(FALSE, TRUE), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

recipe_reference <- function() {
  lapply(seq_len(nrow(recipe_cases)), function(i) {
    m <- do.call(recipe_model, recipe_cases[i, ])
    run(m, ndays = 20L, seed = 1912L)
    recipe_snapshot(m)
  })
}

# Version-2 XDR serialization makes these fingerprints independent of host
# endianness and avoids shipping large simulated histories with the package.
recipe_fingerprint <- function(snapshot) {
  path <- tempfile()
  on.exit(unlink(path))
  writeBin(serialize(snapshot, NULL, version = 2), path)
  unname(tools::md5sum(path))
}
