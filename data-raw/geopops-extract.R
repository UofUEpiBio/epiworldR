# Builds the GeoPops extract shipped in inst/extdata/geopops.
#
# Input: a GeoPops population directory containing
#
#   people_all.csv  net_h.csv  net_w.csv  net_s.csv  net_g.csv
#
# The example used here is Spartanburg County, SC (FIPS 45083): 356,923 people
# and 1,223,555 unique ties across the four layers. That is far too large to
# ship in an R package, so we keep the county's largest census tract plus every
# person directly tied to a resident of it, which leaves the workplaces and
# schools of those residents intact.
#
# Usage:
#   Rscript data-raw/geopops-extract.R <path-to-geopops-population>

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
src  <- if (length(args)) args[[1L]] else stop("Need the GeoPops directory.")
out  <- file.path("inst", "extdata", "geopops")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

layers <- c(h = "net_h", w = "net_w", s = "net_s", g = "net_g")

people <- fread(file.path(src, "people_all.csv"))
nets   <- lapply(layers, function(n)
  fread(file.path(src, paste0(n, ".csv")), select = c("p1", "p2"))
)
edges <- rbindlist(nets, idcol = "layer")

# Seed tract, then one step of closure in the contact network.
tract_sel <- people[tract != 0, .N, by = tract][order(-N)][1, tract]
seed      <- people[tract == tract_sel, uid]
nbrs      <- unique(c(edges[p1 %in% seed, p2], edges[p2 %in% seed, p1]))
keep      <- sort(unique(c(seed, nbrs)))

message(
  "tract ", format(tract_sel, scientific = FALSE),
  ": ", length(seed), " residents -> ", length(keep), " agents"
)

# `sample_index` (an index into the source microdata) and `cbg_geocode`
# (= tract * 10 + block group) are dropped to keep the extract small.
drop_cols <- c("sample_index", "cbg_geocode")

fwrite(
  people[uid %in% keep, setdiff(names(people), drop_cols), with = FALSE],
  file.path(out, "people_all.csv.gz")
)

for (n in names(layers)) {
  e <- nets[[n]][p1 %in% keep & p2 %in% keep]
  e[, edge_weight := 1L]
  fwrite(e, file.path(out, paste0(layers[[n]], ".csv.gz")))
  message(layers[[n]], ": ", nrow(e), " edges")
}
