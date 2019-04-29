# salsa version 0.0.4.900

## changes to R-shiny app

* Cell x-axis can be log-scaled.
* Fit method defaults to "mge" with goodness of fit "KS" for
Frechet-Weibull distribution. The previous "mle" is available
for comparison and testing.

## changes to existing functions

* The fr_wei distribution functions `qfr_wei`, `pfr_wei`, `dfr_wei`
have had the argument `log=FALSE` removed, under guidance of the
`fitdistrplus` package vignette.

## Todo

* evaluate alternative optimization functions, influenced by a blog
post by one of the original authors: https://www.r-bloggers.com/why-optim-is-out-of-date/
who suggests considering the `"optextras"` R package on CRAN.
Major benefit would be substantial speed improvement, and potential
increase in robustness.

# salsa version 0.0.3.900

## minor fixes

* Reference to local debugging files were removed from R-shiny app.

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

