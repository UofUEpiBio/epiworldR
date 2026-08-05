# GeoPops extract

A small extract of the synthetic population from the
[GeoPops Spartanburg measles tutorial](https://github.com/GeoPopsHub/sc_spartanburg_measles/blob/1c741a362dc021b6df082d689f7f7c509f4c23b4/README.md),
used by the vignette "Using external network data: GeoPops synthetic
populations".

The source population is Spartanburg County, SC (FIPS 45083): 356,923 people
and 1,223,555 unique ties. This extract keeps the county's largest census
tract (45083023001) plus every person directly tied to a resident of it, so
that workplaces and schools are not truncated -- 29,534 agents and 90,175
unique ties. See `data-raw/geopops-extract.R` for the exact code.

Because this is a subset, results computed from it are illustrative only, not
estimates for Spartanburg County: agents on the boundary have artificially
truncated neighborhoods.

## Files

| File                | Contents                                                  |
|---------------------|-----------------------------------------------------------|
| `people_all.csv.gz` | One row per person; `uid` is the unique person identifier |
| `net_h.csv.gz`      | Household ties (`p1`, `p2`, `edge_weight`)                |
| `net_w.csv.gz`      | Workplace ties                                            |
| `net_s.csv.gz`      | School ties                                               |
| `net_g.csv.gz`      | Group-quarters ties                                       |

The GeoPops column layout is preserved, except that `sample_index` (an index
into the source microdata) and `cbg_geocode` (`tract * 10 + block group`) were
dropped from `people_all.csv.gz` to keep the files small.

`uid`s in this extract are a sparse subset of the county's, and people who
work or attend school inside the area but live outside it carry no
demographics -- both are handled in the vignette.
