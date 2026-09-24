# Dependencies

Recorded from the environment the 188-loci release was built in. These are the
versions actually installed, not a resolved lockfile — `renv` is not used in
this project, so pinning is by record rather than by restore.

R 4.5.3

| package | version |
|---|---|
| bslib | 0.12.0 |
| circlize | 0.4.18 |
| data.table | 1.18.6.1 |
| dplyr | 1.2.1 |
| DT | 0.34.0 |
| ggplot2 | 4.0.3 |
| openxlsx | 4.2.9 |
| parallel | 4.5.3 |
| pecotmr | 0.8.2 |
| plotly | 4.12.1 |
| readr | 2.2.0 |
| scales | 1.4.0 |
| shiny | 1.14.0 |
| shinycssloaders | 1.1.0 |
| stringr | 1.6.0 |
| tidyverse | 2.0.0 |

`pecotmr` is not on CRAN; see https://github.com/StatFunGen/pecotmr.

One note on `pecotmr`: the locus-level build reads a precomputed variant
correlation table rather than calling `load_LD_matrix()`, which is absent from
the installed version.
