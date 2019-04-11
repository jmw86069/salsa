

#' The Frechet-Weibull Distribution
#'
#' The Frechet-Weibull Distribution
#'
#' Density, distribution function, quantile function,
#' and random generation for the combined Frechet-Weibull
#' distribution with five parameters: `fr_shape`, `fr_scale`,
#' `wei_shape`, `wei_scale`, and `fr_weight`.
#'
#' @family SALSA distribution fit functions
#'
#' @param x vector of quantiles
#' @param fr_shape the Frechet shape
#' @param fr_scale the Frechet scale
#' @param wei_shape the Weibull shape
#' @param wei_scale the Weibull scale
#' @param fr_weight numeric value from 0 to 1, the Frechet weight, where
#'    the corresponding Weibull weight is `(1 - fr_weight)`.
#' @param log logical passed to `actuar::dinvweibull`, if `TRUE`
#'    probabilities/densities are returned as `log(p)`.
#'
#' @export
dfr_wei <- function
(x,
 fr_shape=2,
 fr_scale=100,
 wei_shape=2,
 wei_scale=1000,
 fr_weight=0.5,
 log=FALSE)
{
   ## Purpose is to model two distributions given proper parameters
   ##
   ## fr_weight is the fraction of signal estimated for Frechet,
   ## value between 0 and 1.
   ##
   ## dinvweibull parameters:
   ## shape
   ## scale=1
   ## rate (by default 1/scale)
   ## log=FALSE
   ##
   ## Weibull parameters
   ## shape
   ## scale=1
   ## log=FALSE
   ##
   if (!suppressPackageStartupMessages(require(actuar))) {
      stop("The actuar package is required.");
   }
   if (!suppressPackageStartupMessages(require(jamba))) {
      stop("The jamba package is required.");
   }
   fr_weight <- jamba::noiseFloor(fr_weight, minimum=0, ceiling=1);
   wei_weight <- 1 - fr_weight;
   verbose <- length(getOption("verbose")) > 0 && getOption("verbose");
   if (verbose) {
      printDebug("d_fr_wei(): ",
         "\nfr_shape:", format(fr_shape, digits=2),
         "\nfr_scale:", format(fr_scale, digits=2),
         "\nfr_weight:", format(fr_weight, digits=2),
         "\nwei_shape:", format(wei_shape, digits=2),
         "\nwei_scale:", format(wei_scale, digits=2));
   }

   ## inverse Weibull density (aka Frechet)
   d_fr <- actuar::dinvweibull(x=x,
      shape=fr_shape,
      scale=fr_scale,
      log=log);

   ## Weibull density
   d_wei <- stats::dweibull(x=x,
      shape=wei_shape,
      scale=wei_scale,
      log=log);

   ## combined density
   d_both <- (fr_weight * d_fr + wei_weight * d_wei);
   return(d_both);
}

#' @rdname dfr_wei
#'
#' @export
qfr_wei <- function
(p,
fr_shape=2,
fr_scale=100,
wei_shape=2,
wei_scale=1000,
fr_weight=0.5,
log=FALSE)
{
   if (!suppressPackageStartupMessages(require(actuar))) {
      stop("The actuar package is required.");
   }
   fr_weight <- jamba::noiseFloor(fr_weight, minimum=0, ceiling=1);
   wei_weight <- 1 - fr_weight;
   verbose <- length(getOption("verbose")) > 0 && getOption("verbose");
   if (verbose) {
      printDebug("q_fr_wei(): ",
         "\nfr_shape:", format(fr_shape, digits=2),
         "\nfr_scale:", format(fr_scale, digits=2),
         "\nfr_weight:", format(fr_weight, digits=2),
         "\nwei_shape:", format(wei_shape, digits=2),
         "\nwei_scale:", format(wei_scale, digits=2));
   }

   ## inverse Weibull density (aka Frechet)
   q_fr <- qinvweibull(p,
      shape=fr_shape,
      scale=fr_scale,
      log=log);

   ## Weibull density
   q_wei <- qweibull(p,
      shape=wei_shape,
      scale=wei_scale,
      log=log);

   ## combined quantile
   q_both <- (fr_weight * q_fr + wei_weight * q_wei);
   return(q_both);
}

#' @rdname dfr_wei
#'
#' @export
pfr_wei <- function
(q,
fr_shape=2,
fr_scale=100,
wei_shape=2,
wei_scale=1000,
fr_weight=0.5,
log=FALSE)
{
   if (!suppressPackageStartupMessages(require(actuar))) {
      stop("The actuar package is required.");
   }
   fr_weight <- noiseFloor(fr_weight, minimum=0, ceiling=1);
   wei_weight <- 1 - fr_weight;
   verbose <- length(getOption("verbose")) > 0 && getOption("verbose");
   if (verbose) {
      printDebug("q_fr_wei(): ",
         "\nfr_shape:", format(fr_shape, digits=2),
         "\nfr_scale:", format(fr_scale, digits=2),
         "\nfr_weight:", format(fr_weight, digits=2),
         "\nwei_shape:", format(wei_shape, digits=2),
         "\nwei_scale:", format(wei_scale, digits=2));
   }

   ## inverse Weibull density (aka Frechet)
   p_fr <- pinvweibull(q,
      shape=fr_shape,
      scale=fr_scale,
      log=log);

   ## Weibull density
   p_wei <- pweibull(q,
      shape=wei_shape,
      scale=wei_scale,
      log=log);

   ## combined probability
   p_both <- (fr_weight * p_fr + wei_weight * p_wei);
   return(p_both);
}

#' Parameter bounds for Frechet-Weibull fit
#'
#' Parameter bounds for Frechet-Weibull fit
#'
#' This function defines a `matrix` with default
#' parameter bounds to restrict the range of values
#' allowed during the Frechet-Weibull distribution fit.
#' Each of the five  parameters are given `start`,
#' `lower`, and `upper` values:
#'
#' * `fr_weight` the weight of Frechet relative to Weibull,
#'    on a scale of 0 to 1.
#' * `fr_shape` the shape value for the Frechet distribution
#' * `fr_scale` the scale value for the Frechet distribution
#' * `wei_shape` the shape value for the Weibull distribution
#' * `wei_scale` the scale value for the Weibull distribution
#'
#' The output `matrix` can be modified then passed to
#' `do_salsa_steps()` in order to customize the model fit.
#'
#' @return `matrix` of five parameter upper, lower, and start
#'    values to be used by `fitdist_fr_wei()`.
#'
#' @family SALSA distribution fit functions
#'
#' @param param_fr_wei `data.frame` containing one or more
#'    colnames `"fr_shape"`, `"fr_scale"`, `"wei_shape"`,
#'    `"wei_scale"`, `"fr_weight"`, and one or more rownames
#'    `"start"`, `"lower"`, `"upper"`. Any non-NA values
#'    replace the corresponding default values.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' # default matrix of parameters
#' param_fr_wei <- params_fr_wei();
#' param_fr_wei;
#'
#' # modify some values
#' custom_param_fr_wei_m <- rbind(
#'    start=c(fr_shape=50, wei_shape=40),
#'    lower=c(fr_shape=NA, wei_shape=400));
#' custom_param_fr_wei_m;
#' custom_param_fr_wei <- params_fr_wei(custom_param_fr_wei_m);
#' custom_param_fr_wei;
#'
#' @export
params_fr_wei <- function
(param_fr_wei=NULL,
   ...)
{
   ## Define a matrix with start, low, and high parameter values
   ## for Frechet-Weibull distribution fit parameters
   lowerFrWei <- c(fr_weight=0.01,
      fr_shape=1.5,
      fr_scale=1,
      wei_shape=1.5,
      wei_scale=1);
   upperFrWei <- c(fr_weight=0.99,
      fr_shape=1000,
      fr_scale=10000,
      wei_shape=1000,
      wei_scale=10000);
   startFrWei <- c(fr_weight=0.1,
      fr_shape=2,
      fr_scale=179,
      wei_shape=2,
      wei_scale=1700);
   param_fr_wei1 <- rbind(start=startFrWei,
      lower=lowerFrWei,
      upper=upperFrWei);
   if (length(param_fr_wei) > 0) {
      # replace any non-NA entries from supplied colnames and rownames
      colv <- intersect(colnames(param_fr_wei1), colnames(param_fr_wei));
      rowv <- intersect(rownames(param_fr_wei1), rownames(param_fr_wei));
      if (length(colv) > 0 && length(rowv) > 0) {
         for (coli in colv) {
            non_na <- !is.na(param_fr_wei[rowv,colv]);
            if (any(non_na)) {
               param_fr_wei1[rowv,colv][non_na] <- param_fr_wei[non_na];
            }
         }
      }
   }
   return(param_fr_wei1);
}

#' Fit Frechet-Weibull distribution to non-censored data
#'
#' Fit Frechet-Weibull distribution to non-censored data
#'
#' This function is a simple wrapper around `fitdistrplus::fitdist()`
#' for the five-parameter Frechet-Weibull distribution.
#' Most customizations are performed by the corresponding
#' support function, for example when `method="mle"` the
#' distribution fit function is `fitdistrplus::mledist()`,
#' which then calls `stats::optim()`. Arguments are
#' passed to these functions using `...` where applicable.
#'
#' TODO: This function currently expects to use `method="mle"`
#' which utilizes the upper,lower,start parameters for each
#' of the five parameters. Other methods may have different
#' requirements which are not formally handled. In future
#' we will evaluate and support additional methods, modifying
#' the call to `fitdistrplus::fitdist()` as needed.
#'
#' @return object of class `"fitdist"` described in
#'     `fitdistrplus::fitdist()`.
#'
#' @family SALSA distribution fit functions
#'
#' @param x numeric vector
#' @param param_fr_wei `data.frame` of parameter bounds, output
#'    from `params_fr_wei()`, which defines the start, upper, and
#'    lower bounds for each of the five parameters.
#' @param method fitting method, argument passed to
#'    `fitdistrplus::fitdist()`.
#' @param optim.method character `method` argument passed to
#'    `fitdistrplus::mledist` or corresponding function based
#'    upon the `method` argument, but ultimately passed to
#'    `stats::optim()`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `fitdistrplus::fitdist()`.
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#' x <- oz2_numi_per_cell$count[oz2_numi_per_cell$count >= 16];
#' fit_fr_wei <- fitdist_fr_wei(x);
#' # print coefficients
#' print(coef(fit_fr_wei));
#'
#' @export
fitdist_fr_wei <- function
(x,
 param_fr_wei=params_fr_wei(),
 method="mle",
 optim.method="Nelder-Mead",
 log=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide a simple wrapper for fitdist()
   if (verbose) {
      jamba::printDebug("fitdist_fr_wei(): ",
         "method:", method);
      jamba::printDebug("fitdist_fr_wei(): ",
         "optim.method:", optim.method);
      jamba::printDebug("fitdist_fr_wei(): ",
         "param_fr_wei:");
      print(param_fr_wei);
   }
   fit_fr_wei <- tryCatch({
      fitdistrplus::fitdist(
         data=x,
         distr="fr_wei",
         method=method,
         start=as.list(param_fr_wei["start",]),
         lower=param_fr_wei["lower",],
         upper=param_fr_wei["upper",],
         optim.method=optim.method,
         ...);
   }, error=function(e){
      if (verbose) {
         jamba::printDebug("fitdist_fr_wei(): ",
            "Error:");
         print(e);
      }
      NULL;
   });
   return(fit_fr_wei);
}

#' Fit Frechet distribution to non-censored data
#'
#' Fit Frechet distribution to non-censored data
#'
#' This function is a simple wrapper around `fitdistrplus::fitdist()`
#' for the Frechet distribution, using `"actuar::dinvweibull"`.
#' Most customizations are performed by the corresponding
#' support function, for example when `method="mme"` the
#' distribution fit function is `fitdistrplus::mmedist()`,
#' which then calls `stats::optim()`. Arguments are
#' passed to these functions using `...` where applicable.
#'
#' @return object of class `"fitdist"` described in
#'     `fitdistrplus::fitdist()`. If any error occurs during
#'     the fit, a `NULL` is returned.
#'
#' @family SALSA distribution fit functions
#'
#' @param x numeric vector
#' @param param_fr_wei `data.frame` of parameter bounds, output
#'    from `params_fr_wei()`, which defines the start, upper, and
#'    lower bounds for each of the five parameters.
#' @param method fitting method, argument passed to
#'    `fitdistrplus::fitdist()`, by default `"mme"` moment
#'    matching estimation.
#' @param min_shape numeric value restricting the minimum shape
#'    otherwise output is set to `NULL`.
#' @param order,central,absolute,na.rm arguments passed to
#'    `moments::moment()`. These values do not typically need to
#'    be customized.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `fitdistrplus::fitdist()`.
#'
#' @examples
#' data(oz2_numi_per_cell);
#' x <- oz2_numi_per_cell$count[oz2_numi_per_cell$count >= 16];
#' fit_fr <- fitdist_fr(x);
#' # print coefficients
#' print(coef(fit_fr));
#'
#' @export
fitdist_fr <- function
(x,
 method="mme",
 min_shape=1,
 order=c(1,2),
 central=FALSE,
 absolute=FALSE,
 na.rm=FALSE,
 verbose=FALSE,
 ...)
{
   ## default moment function
   get_memp2 <- function(
      central=FALSE,
      absolute=FALSE,
      na.rm=FALSE) {
      memp2 <- function(x, order) {
         moments::moment(x,
            order,
            central,
            absolute,
            na.rm);
      }
      return(memp2);
   }
   memp2 <- get_memp2(central, absolute, na.rm);

   ## fitdist but return NULL if error occurs
   fit_fr <- tryCatch({
      fitdistrplus::fitdist(x,
         distr="invweibull",
         method=method,
         order=order,
         memp=memp2);
   }, error=function(e) {
      if (verbose) {
         printDebug("fitdist_fr(): ",
            "Error:");
         print(e);
      }
      NULL;
   });
   if (length(fit_fr) > 0 &&
         !is.na(fit_fr$estimate["shape"]) &&
         length(min_shape) > 0 &&
         fit_fr$estimate["shape"] >= min_shape) {
      fit_fr;
   } else {
      NULL;
   }
}

#' Calculate lower bound from Frechet-Weibull fit
#'
#' Calculate lower bound from Frechet-Weibull fit
#'
#' This function takes the output from `fitdist_fr_wei()` and
#' calculates the associated predicted lower bound for
#' number of UMI to accept for downstream analysis of
#' scRNA-seq data.
#'
#' The equation:
#'
#' fr_scale * (wei_scale / fr_scale) ^ (fr_weight / (fr_shape * wei_shape))
#'
#' @family SALSA support functions
#'
#' @return numeric value representing the lowest number of UMI to
#' accept for downstream analysis of scRNA-seq data.
#'
#' @param fit_fr_wei `fitdist` object output from `fitdist_fr_wei()`,
#'    or numeric vector with names `"fr_weight"`, `"fr_shape"`,
#'    `"fr_scale"`, `"wei_shape"`, `"wei_scale"`.
#' @param ... additional arguments are ignored.
#'
#' @export
get_lower_bound <- function
(fit_fr_wei,
 ...)
{
   ## Purpose is to calculate the recommended lower bound
   ## based upon a Frechet-Weibull fit result.
   if ("numeric" %in% class(fit_fr_wei)) {
      coefx <- c("fr_weight", "fr_shape", "fr_scale",
         "wei_shape", "wei_scale");
      if (!all(coefx %in% names(fit_fr_wei))) {
         stop(paste0("Numeric input must include all names: ",
            paste(coefx, collapse=",")));
      }
      coef_fr_wei <- fit_fr_wei;
   } else {
      coef_fr_wei <- coef(fit_fr_wei);
   }
   lower_bound <- coef_fr_wei["fr_scale"] *
      ( coef_fr_wei["wei_scale"] / coef_fr_wei["fr_scale"] ) ^
      ( coef_fr_wei["fr_weight"] /
            ( coef_fr_wei["fr_shape"] * coef_fr_wei["wei_shape"] ) );
   return(unname(lower_bound));
}

#' Calculate upper bound from Frechet fit
#'
#' Calculate upper bound from Frechet fit
#'
#' This function takes the output from `fitdist_fr()` and
#' calculates the associated predicted upper bound for
#' number of UMI to accept for downstream analysis of
#' scRNA-seq data.
#'
#' The equation:
#'
#' scale + 2 * shape * scale * (base::gamma(1 - (1 / shape)) - 1)
#'
#' @return numeric value representing the highest number of UMI to
#' accept for downstream analysis of scRNA-seq data.
#'
#' @family SALSA support functions
#'
#' @param fit_fr `fitdist` object output from `fitdist_fr()`,
#'    or numeric vector with names `"shape"`, `"scale"`.
#' @param ... additional arguments are ignored.
#'
#' @export
get_upper_bound <- function
(fit_fr,
   ...)
{
   ## Purpose is to calculate the recommended lower bound
   ## based upon a Frechet fit result.
   if ("numeric" %in% class(fit_fr)) {
      coefx <- c("shape", "scale");
      if (!all(coefx %in% names(fit_fr))) {
         stop(paste0("Numeric input must include all names: ",
            paste(coefx, collapse=",")));
      }
      coef_fr <- fit_fr;
   } else {
      coef_fr <- coef(fit_fr);
   }
   upper_bound <- (coef_fr["scale"]) +
      2 * coef_fr["shape"] * coef_fr["scale"] *
      (base::gamma(1 - (1 / coef_fr["shape"])) - 1);

}


#' Get step parameters for SALSA
#'
#' Get step parameters for SALSA
#'
#' This function takes a vector of counts and determines the
#' appropriate step size to use when iterating the count
#' threshold used by SALSA.
#'
#' @return `list` containing `"n_start"`, `"step_size"`,
#'    `"max_step"`, and when `include_vector=TRUE` it
#'    includes `"n_vector"` and `"count_vector"`.
#'
#' @family SALSA support functions
#'
#' @param x numeric vector of counts
#' @param ... additional arguments are ignored
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#'
#' usecounts <- sort(oz2_numi_per_cell$count);
#' get_salsa_steps(usecounts);
#'
#' # optionally return the vector of thresholds to use
#' get_salsa_steps(usecounts, include_vector=TRUE);
#'
#' @export
get_salsa_steps <- function
(x,
 include_vector=FALSE,
 ...)
{
   #
   minUMI <- min(x, na.rm=TRUE);
   maxUMI <- max(x, na.rm=TRUE);
   n_start <- max(c(floor(sqrt(minUMI)), 1));
   n_bound <- floor(sqrt(maxUMI));
   step_size <- ceiling( log10 ( n_bound ) );
   max_step <- floor(maxUMI^(1/step_size));

   ret_vals <- list(n_start=n_start,
      step_size=step_size,
      max_step=max_step);
   if (include_vector) {
      n_vector <- seq(from=n_start,
         to=max_step,
         by=1);
      count_vector <- (n_vector ^ step_size) - 1;
      ## I also saw this variation on the equation, seems incorrect
      #count_vector <- (n_vector-1) ^ step_size;
      ret_vals$n_vector <- n_vector;
      ret_vals$count_vector <- count_vector;
   }
   return(ret_vals);
}

#' Perform SALSA steps for threshold detection
#'
#' Perform SALSA steps for threshold detection
#'
#' This function is a wrapper around `fitdist_fr()` and
#' `fitdist_fr_wei()`, which iterates through a wide range
#' of possible thresholds to determine the fit parameters,
#' and associated lower and upper bounds. The results
#' are intended to be plotted to determine appropriate
#' thresholds to use when calculating the lower and upper
#' bounds for barcodes and genes in a single cell RNA-seq
#' dataset.
#'
#' @param x numeric vector of counts, either the number of
#'    UMI per cell, or the number of UMI per gene.
#' @param n_vector,n_start,max_step,step_size,count_vector
#'    arguments passed to `get_salsa_steps()` which returns
#'    `count_vector`. If `count_vector` is supplied, it is
#'    used without modification.
#' @param dists character vector determining which distribution
#'    fit functions to calculate, `"frechet"` fits the Frechet
#'    distribution, `"frechet-weibull"` fits the five-parameter
#'    combined Frechet and Weibull distributions. Typically the
#'    Frechet parameters are used to define the upper bound,
#'    and the Frechet-Weibull parameters are used to define the
#'    lower bound.
#' @param cache_fr,cache_fr_wei list objects output from
#'    `memoise::cache_filesystem()` used to store cached
#'    distribution fit results. When either are NULL, the
#'    memoise cache steps are disabled, and the functions are
#'    called directly.
#' @param param_fr_wei `data.frame` output from `params_fr_wei()`
#'    which defines the parameter upper and lower limits, and
#'    start values for each parameter. If `NULL` then the
#'    default values from `params_fr_wei()` are used.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @return `list` with one element for each value in `count_vector`,
#'    where each list element contains a list with one entry for
#'    each value in argument `dists` containing the fit parameters
#'    for each selected distribution, as well as an entry `"min_count"`
#'    which contains the minimum counts to use in each fit. When
#'    `dists` contains `"frechet-weibull"` each list includes
#'    `"lower_bound"`. When `dists` contains `"frechet"` each list
#'    includes `"upper_bound"`. The output is intended to be passed
#'    to `get_salsa_table()`.
#'
#' @family SALSA core functions
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#' x <- oz2_numi_per_cell$count[oz2_numi_per_cell$count >= 16];
#' x_salsa <- do_salsa_steps(x,
#'    count_vector=c(16,32,128),
#'    cache_fr=NULL,
#'    cache_fr_wei=NULL);
#' x_df <- get_salsa_table(x_salsa);
#' x_df;
#'
#' @export
do_salsa_steps <- function
(x,
 n_vector=NULL,
 n_start=NULL,
 max_step=NULL,
 step_size=NULL,
 count_vector=NULL,
 dists=c("frechet", "frechet-weibull"),
 cache_fr=cache_filesystem("./cache_fr"),
 cache_fr_wei=cache_filesystem("./cache_fr_wei"),
 param_fr_wei=NULL,
 verbose=FALSE,
 ...)
{
   ## Set up memoise for caching
   if (length(cache_fr) == 0) {
      fitdist_fr_m <- fitdist_fr;
   } else {
      fitdist_fr_m <- memoise::memoise(fitdist_fr,
         cache=cache_fr);
   }
   if (length(cache_fr) == 0) {
      fitdist_fr_wei_m <- fitdist_fr_wei;
   } else {
      fitdist_fr_wei_m <- memoise::memoise(fitdist_fr_wei,
         cache=cache_fr);
   }

   # verify count_vector by several possible methods
   if (length(count_vector) == 0) {
      ## count_vector not supplied so it must be created
      if (length(n_vector) == 0) {
         if (length(n_start) == 0 ||
               length(n_bound) == 0 ||
               length(step_size) == 0) {
            step_list <- get_salsa_steps(x,
               include_vector=TRUE);
         }
         if (length(n_start) == 0) {
            n_start <- step_list$n_start;
         }
         if (length(max_step) == 0) {
            max_step <- step_list$max_step;
         }
         if (length(step_size) == 0) {
            step_size <- step_list$step_size;
         }
         n_vector <- seq(from=n_start,
            to=max_step,
            by=1);
      }
      if (verbose) {
         jamba::printDebug("do_salsa_steps(): ",
            "n_start:", n_start);
         jamba::printDebug("do_salsa_steps(): ",
            "step_size:", step_size);
         jamba::printDebug("do_salsa_steps(): ",
            "max_step:", max_step);
         jamba::printDebug("do_salsa_steps(): ",
            "n_vector:", n_vector);
      }
      count_vector <- (n_vector ^ step_size) - 1;
   }
   if (verbose) {
      jamba::printDebug("do_salsa_steps(): ",
         "count_vector:",
         count_vector);
   }

   ## Define Frechet-Weibull starting parameters
   if (length(param_fr_wei) == 0) {
      if (verbose) {
         jamba::printDebug("do_salsa_steps(): ",
            "using defaults for ",
            "param_fr_wei");
      }
      param_fr_wei <- params_fr_wei();
   }

   ## Iterate each n_vector value:
   ## - use memoise(fitdist_fr_wei()) to cache Frechet-Weibull fit results
   ## - use memoise(fitdist_fr()) to cache Frechet fit results
   fit_list <- lapply(jamba::nameVector(count_vector), function(icount){
      if (verbose) {
         jamba::printDebug("do_salsa_steps(): ",
            "icount:", icount);
      }
      ret_vals <- list();
      ret_vals$min_count <- icount;
      usecounts <- x[x >= icount];
      if (jamba::igrepHas("^frechet$", dists)) {
         ## Fit Frechet
         if (verbose) {
            jamba::printDebug("do_salsa_steps(): ",
               "fitdist_fr_m()");
         }
         fit_fr <- tryCatch({
            fitdist_fr_m(usecounts,
               ...);
         }, error=function(e){
            if (verbose) {
               jamba::printDebug("do_salsa_steps(): ",
                  "fitdist_fr_m() error:");
               print(e);
            }
            NULL;
         });
         #fit_fr_v <- as.list(fitFL$estimate);
         ret_vals$fit_fr <- fit_fr;
         # get upper bound
         if (length(fit_fr) > 0) {
            upper_bound <- get_upper_bound(fit_fr);
         } else {
            upper_bound <- NA;
         }
         ret_vals$upper_bound <- upper_bound;
      }
      if (jamba::igrepHas("frechet-weibull", dists)) {
         if (length(usecounts) < 6) {
            ret_vals$fit_fr_wei <- NULL;
            ret_vals$lower_bound <- NA;
         } else {
            ## Fit Frechet
            if (verbose) {
               jamba::printDebug("do_salsa_steps(): ",
                  "fitdist_fr_wei_m()");
            }
            fit_fr_wei <- tryCatch({
               fitdist_fr_wei_m(x=usecounts,
                  param_fr_wei=param_fr_wei);
            }, error=function(e){
               if (verbose) {
                  jamba::printDebug("do_salsa_steps(): ",
                     "fitdist_fr_wei_m() error:");
                  print(e);
               }
               NULL;
            });
            ret_vals$fit_fr_wei <- fit_fr_wei;
            ## get lower bound
            if (length(fit_fr_wei) > 0) {
               lower_bound <- get_lower_bound(fit_fr_wei);
            } else {
               lower_bound <- NA;
            }
            ret_vals$lower_bound <- lower_bound;
         }
      }
      ret_vals;
   });
   return(fit_list);
}

#' Get SALSA results as a data.frame
#'
#' Get SALSA results as a data.frame
#'
#' This function converts the list of lists output from
#' `do_salsa_steps()` into a "tidy" data.frame intended
#' to be used for visualization.
#'
#' @return `data.frame` containing one row per minimum count
#'    threshold, with columns containing the Frechet and
#'    Frechet-Weibull fit parameters, or NA values when there
#'    are no parameters for the given minimum count threshold.
#'
#' @family SALSA core functions
#'
#' @param fit_list `list` output from `do_salsa_steps()`,
#'    where each list element contains `min_count`,
#'    and one or more of `fit_fr` and `fit_fr_wei`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#' x <- oz2_numi_per_cell$count[oz2_numi_per_cell$count >= 16];
#' x_salsa <- do_salsa_steps(x,
#'    count_vector=c(16,32,128),
#'    cache_fr=NULL,
#'    cache_fr_wei=NULL);
#' x_df <- get_salsa_table(x_salsa);
#' x_df;
#'
#' @export
get_salsa_table <- function
(fit_list,
 verbose=FALSE,
 ...)
{
   ##
   if (length(fit_list) == 0) {
      return(NULL);
   }
   fit_list_df <- rbindList(lapply(nameVectorN(fit_list), function(icount){
      n <- fit_list[[icount]]$min_count;
      if (verbose) {
         printDebug("get_salsa_table(): ",
            "icount:", icount);
      }
      if ("fit_fr" %in% names(fit_list[[icount]])) {
         fr_coef <- as.list(coef(fit_list[[icount]]$fit_fr));
      } else {
         fr_coef <- list(shape=NA,
            scale=NA);
      }
      if ("fit_fr_wei" %in% names(fit_list[[icount]])) {
         fr_wei_coef <- as.list(coef(fit_list[[icount]]$fit_fr_wei));
      } else {
         fr_wei_coef <- list(fr_weight=NA,
            fr_shape=NA,
            fr_scale=NA,
            wei_shape=NA,
            wei_scale=NA);
      }
      data.frame(count=n,
         shape=fr_coef$shape,
         scale=fr_coef$scale,
         fr_weight=fr_wei_coef$fr_weight,
         fr_shape=fr_wei_coef$fr_shape,
         fr_scale=fr_wei_coef$fr_scale,
         wei_shape=fr_wei_coef$wei_shape,
         wei_scale=fr_wei_coef$wei_scale,
         lower_bound=fit_list[[icount]]$lower_bound,
         upper_bound=fit_list[[icount]]$upper_bound
      );
   }));
   fit_list_df;
}
