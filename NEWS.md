# salsa version 0.0.2.900

## additions

* First prototype of R-shiny app for the nUMI per cell and gene
threshold workflow.

# salsa version 0.0.1.900

## additions

* Added examples for several functions.
* Two example datasets are included to help illustrate the workflow,
`oz2_numi_per_cell` and `oz2_numi_per_gene`.
* `get_salsa_steps()`, `do_salsa_steps()` and `get_salsa_table()` are
the core workflow for determining the upper and lower bound for
scRNA-seq barcode and gene values to use in downstream analyses.

## changes to existing functions

* `do_salsa_steps()` uses memoise to cache result for distribution fit
functions `fitdist_fr()` and `fitdist_fr_wei()`.
* `do_salsa_steps()` accepts custom parameter boundaries for the
Frechet-Weibull fit, in `data.frame` format as returned by
`params_fr_wei()`. The recommended workflow is to generate a default
`data.frame` from `params_fr_wei()`, modify its values, then
pass it as an argument to `do_salsa_steps()`.
* `params_fr_wei()` now accepts partial custom input, and replaces
non-NA values for only the rows and columns supplied, making it easier
to supply only certain custom parameters.

# salsa version 0.0.0.900

* Initial package creation.

