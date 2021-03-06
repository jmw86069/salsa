% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/salsa-plot.R
\name{plot_fr_wei}
\alias{plot_fr_wei}
\title{Plot Frechet-Weibull fit}
\usage{
plot_fr_wei(fit, addx = NULL, ylim = NULL, xlim = NULL,
  scale_y = FALSE, log = "x", type = c("gg", "base"),
  highlight_range = NULL, do_facet = TRUE, ...)
}
\arguments{
\item{fit}{\code{fitdist} object as output from \code{fitdist_fr_wei()}.
Optionally \code{fit} can be a named vector with \code{"fr_weight"},
\code{"fr_shape"}, \code{"fr_scale"}, \code{"wei_shape"}, \code{"wei_scale"}.
In that case, use \code{addx} to supply numeric range of data
to plot.}

\item{addx}{optional numeric vector of x values to add to
the x-axis range.}

\item{ylim, xlim}{optional numeric vector describing the y-axis
or x-axis limits, respectively.}

\item{scale_y}{logical indicating whether to scale the densities
to common maximum height of y=1.}

\item{log}{character vector indicating the axis or axes
log-transform, with any presence of \code{"x"} and \code{"y"}.}

\item{type}{character string indicating the type of plot,
where \code{"gg"} uses ggplot2, and \code{"base"} uses base R
graphics.}

\item{highlight_range}{optional numeric range to highlight,
intended to represent the range of counts selected for
downstream analysis.}

\item{do_facet}{logical indicating whether to create a
\code{ggplot2::facet_wrap()} output, which puts a colored
title bar above the plot panel.}

\item{...}{additional arguments are ignored.}
}
\description{
Plot Frechet-Weibull fit
}
\details{
This function plots the two distributions, Frechet and Weibull,
using the fit parameters provided.
}
\examples{
library(salsa);
data(oz2_numi_per_cell);
param_fr_wei <- params_fr_wei();
usecounts <- sort(oz2_numi_per_cell$count);
usecounts <- usecounts[usecounts >= 16];
fit1 <- fitdist_fr_wei(x=usecounts, param_fr_wei=param_fr_wei);
plot_fr_wei(fit1);
plot_fr_wei(fit1, ylim=c(0,0.005));
plot_fr_wei(fit1, scale_y=TRUE) + ggtitle("densities scaled to max height")
plot_fr_wei(fit1, scale_y=TRUE, highlight_range=c(100, 2200))

}
\seealso{
Other SALSA visualization functions: \code{\link{rank_count_plot}}
}
\concept{SALSA visualization functions}
