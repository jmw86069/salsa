
#' Fit Frechet-Weibull distributions
#'
#' Fit Frechet-Weibull distributions
#'
#' @export
fit_fr_wei <- function
(...)
{
   ## Purpose is to fit the five-parameter distribution
   ## combining Weibull and Frechet distributions.
}

#' The Frechet-Weibull Distribution
#'
#' The Frechet-Weibull Distribution
#'
#' Density, distribution function, quantile function,
#' and random generation for the combined Frechet-Weibull
#' distribution with five parameters: `fr_shape`, `fr_scale`,
#' `wei_shape`, `wei_scale`, and `fr_weight`.
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
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
dfr_wei <- function
(x,
 fr_shape=2,
 fr_scale=100,
 wei_shape=2,
 wei_scale=1000,
 fr_weight=0.5,
 log=FALSE,
 verbose=FALSE,
 ...)
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
log=FALSE,
verbose=FALSE,
...)
{
   if (!suppressPackageStartupMessages(require(actuar))) {
      stop("The actuar package is required.");
   }
   fr_weight <- jamba::noiseFloor(fr_weight, minimum=0, ceiling=1);
   wei_weight <- 1 - fr_weight;
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
log=FALSE,
verbose=FALSE,
...)
{
   if (!suppressPackageStartupMessages(require(actuar))) {
      stop("The actuar package is required.");
   }
   fr_weight <- noiseFloor(fr_weight, minimum=0, ceiling=1);
   wei_weight <- 1 - fr_weight;
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
#' This function defines a data.frame of parameter bounds,
#' restricting the range of values allowed during the
#' Frechet-Weibull distribution fit. Each of the five
#' parameters are given `start`, `lower`, and `upper`
#' values:
#'
#' * `fr_weight` the weight of Frechet relative to Weibull,
#'    on a scale of 0 to 1.
#' * `fr_shape` the shape value for the Frechet distribution
#' * `fr_scale` the scale value for the Frechet distribution
#' * `wei_shape` the shape value for the Weibull distribution
#' * `wei_scale` the scale value for the Weibull distribution
#'
#' @param param_fr_wei `data.frame` containing colnames
#'    `"fr_shape"`, `"fr_scale"`, `"wei_shape"`, `"wei_scale"`,
#'    and `"fr_weight"`, and rownames `"start"`, `"lower"`,
#'    `"upper"`.
#' @param ... additional arguments are ignored.
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
      colx <- c("fr_weight","fr_shape", "fr_scale", "wei_shape", "wei_scale");
      rowx <- c("start", "lower", "upper");
      param_fr_wei <- param_fr_wei[rowx,colx];
      if (any(!is.na(param_fr_wei))) {
         param_fr_wei1[!is.na(param_fr_wei)] <- param_fr_wei[!is.na(param_fr_wei)];
      }
   }
   param_fr_wei <- param_fr_wei1;
   return(param_fr_wei);
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
#' @return object of class `"fitdist"` described in
#'     `fitdistrplus::fitdist()`.
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
#' @export
fitdist_fr_wei <- function
(x,
 param_fr_wei,
 method="mle",
 optim.method="Nelder-Mead",
 log=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide a simple wrapper for fitdist()
   if (verbose) {
      printDebug("fitdist_fr_wei(): ",
         "method:", method);
      printDebug("fitdist_fr_wei(): ",
         "optim.method:", optim.method);
      printDebug("fitdist_fr_wei(): ",
         "param_fr_wei:");
      print(param_fr_wei);
   }
   fit_fr_wei <- fitdistrplus::fitdist(
      data=x,
      distr="fr_wei",
      method=method,
      start=as.list(param_fr_wei["start",]),
      lower=param_fr_wei["lower",],
      upper=param_fr_wei["upper",],
      optim.method=optim.method,
      ...);
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
#'    `moments::moment()`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `fitdistrplus::fitdist()`.
#'
#' @examples
#' data(oz2_numi_per_cell);
#' usecounts <- oz2_numi_per_cell$count[oz2_numi_per_cell$count >= 16];
#' fit_fr <- fitdist_fr(x=usecounts);
#' coef(fit_fr);
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
      printDebug("Error:");
      print(e);
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


fit_weibull_test <- function
(x)
{
   xbar <- mean(x)
   varx <- var(x)
   f <- function(b){
      gamma(1+2/b) /  gamma(1+1/b)^2 - 1 - varx/xbar^2
   }
   bhat <- uniroot(f,c(0.02,50))$root
   ahat <- xbar/gamma(1+1/bhat);
   return(c(ahat,bhat))
}

fitdist_weibull_test <- function
(x)
{
   fitdistrplus::fitdist(x,
      distr="weibull",
      method="mle");
}

#' Get step parameters for SALSA
#'
#' Get step parameters for SALSA
#'
#' This function takes a vector of counts and determines the
#' appropriate step size to use when iterating the count
#' threshold used by SALSA.
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
   step_size;
   ret_vals <- list(n_start=n_start,
      step_size=step_size,
      n_bound=n_bound);
   if (include_vector) {
      n_vector <- seq(from=n_start,
         to=n_bound,
         by=step_size);
      count_vector <- (n_vector^2) - 1;
      ret_vals$n_vector <- n_vector;
      ret_vals$count_vector <- count_vector;
   }
   return(ret_vals);
}

#' Perform SALSA steps for threshold detection
#'
#' Perform SALSA steps for threshold detection
#'
#' @export
do_salsa_steps <- function
(x,
 n_vector=NULL,
 n_start=NULL,
 n_bound=NULL,
 step_size=NULL,
 ...)
{
   #
   if (length(n_vector) == 0) {
      if (length(n_start) == 0 ||
            length(n_bound) == 0 ||
            length(step_size) == 0) {
         step_list <- get_salsa_steps(x, include_vector=TRUE);
      }
      if (length(n_start) == 0) {
         n_start <- step_list$n_start;
      }
      if (length(n_bound) == 0) {
         n_bound <- step_list$n_bound;
      }
      if (length(step_size) == 0) {
         step_size <- step_list$step_size;
      }
      n_vector <- seq(from=n_start,
         to=n_bound,
         by=step_size);
   }
   n_vector;
   ## Iterate each n_vector value:
   ## - use memoise(fitdist_fr_wei()) to cache Frechet-Weibull fit results
   ## - use memoise(fitdist_fr()) to cache Frechet fit results
}
