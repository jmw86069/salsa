
#' Rank-Count plot
#'
#' Rank-Count plot
#'
#' This function creates a classic "knee plot" showing
#' `"Rank"` on the x-axis, and `"Count"` on the y-axis,
#' typically with log-scaled x- and y-axes.
#'
#' @param x numeric or integer vector of counts
#' @param log character vector indicating the axis or axes
#'    log-transform, with any presence of `"x"` and `"y"`.
#' @param flip_axes logical indicating whether to flip the
#'    x- and y-axis coordinates.
#' @param highlight_range optional numeric range to highlight,
#'    intended to represent the range of counts selected for
#'    downstream analysis.
#' @param ... additional arguments are ignored
#'
#' @return `ggplot` object sufficient to plot the Rank-Count
#'    data, including the optional highlighted points based
#'    upon the `highlight_range` values.
#'
#' @family SALSA visualization functions
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#'
#' usecounts <- sort(oz2_numi_per_cell$count);
#'
#' rank_count_plot(usecounts,
#'    highlight_range=c(100,2500)) +
#'    ggtitle("All counts unfiltered");
#'
#' # after filtering, it may help add y=1 to the axis range
#' usecounts <- usecounts[usecounts >= 16];
#' rank_count_plot(usecounts,
#'    highlight_range=c(100,2500)) +
#'    expand_limits(y=c(1,NA)) +
#'    ggtitle("Filtered for counts >= 16");
#'
#' @export
rank_count_plot <- function
(x,
 log=c("x", "y"),
 flip_axes=FALSE,
 highlight_range=NULL,
 ...)
{
   ## Purpose is to take CellData which contains rank and counts,
   ## and produce a log-scaled ggplot object.
   df <- data.frame(Count=sort(x));
   df$Rank <- rank(ties.method="min",
      -df$Count);

   ## make unique to remove duplicated rank values
   df <- unique(df);

   ## Optionally highlight a range of points
   if (length(highlight_range) > 0) {
      highlight_flag <- df$Count >= min(highlight_range) & df$Count <= max(highlight_range);
      df$highlight <- ifelse(highlight_flag,
         "selected", "discarded");
   } else {
      df$highlight <- "discarded";
   }

   if (flip_axes) {
      gg <- ggplot(data=df,
         aes(x=Count,
            y=Rank));
   } else {
      gg <- ggplot(data=df,
         aes(x=Rank,
            y=Count));
   }
   gg <- gg +
      geom_line(show.legend=FALSE, alpha=0.5) +
      geom_point(show.legend=FALSE, aes(color=highlight)) +
      #facet_grid(type~name, scales="free_y")+
      colorjam::theme_jam() +
      colorjam::scale_color_jam();
   if (length(highlight_range) > 0) {
      gg <- gg +
         ggforce::geom_mark_rect(#data=subset(df, highlight %in% "selected"),
            aes(filter=highlight=="selected",
               label=highlight),
            expand=unit(0.01, "mm"),
            radius=unit(0.01, "mm"));
   }
   ## Add x-axis log scale
   if (igrepHas("x", log)) {
      gg <- gg +
         scale_x_log10(
            labels=function(x) {
               jamba::asSize(x,
                  byteSize=1000,
                  unitType="",
                  sep="");
            }) +
         annotation_logticks(base=10);
   }
   if (igrepHas("y", log)) {
      gg <- gg +
         scale_y_log10(
            labels=function(x) {
               jamba::asSize(x,
                  byteSize=1000,
                  unitType="",
                  sep="");
            }) +
         annotation_logticks(base=10);
   }
   return(gg);
}

#' Plot Frechet-Weibull fit
#'
#' Plot Frechet-Weibull fit
#'
#' This function plots the two distributions, Frechet and Weibull,
#' using the fit parameters provided.
#'
#' @param fit `fitdist` object as output from `fitdist_fr_wei()`.
#' @param addx optional numeric vector of x values to add to
#'    the x-axis range.
#' @param ylim,xlim optional numeric vector describing the y-axis
#'    or x-axis limits, respectively.
#' @param scale_y logical indicating whether to scale the densities
#'    to common maximum height of y=1.
#' @param log character vector indicating the axis or axes
#'    log-transform, with any presence of `"x"` and `"y"`.
#' @param type character string indicating the type of plot,
#'    where `"gg"` uses ggplot2, and `"base"` uses base R
#'    graphics.
#' @param highlight_range optional numeric range to highlight,
#'    intended to represent the range of counts selected for
#'    downstream analysis.
#' @param do_facet logical indicating whether to create a
#'    `ggplot2::facet_wrap()` output, which puts a colored
#'    title bar above the plot panel.
#' @param ... additional arguments are ignored.
#'
#' @family SALSA visualization functions
#'
#' @examples
#' library(salsa);
#' data(oz2_numi_per_cell);
#' param_fr_wei <- params_fr_wei();
#' usecounts <- sort(oz2_numi_per_cell$count);
#' usecounts <- usecounts[usecounts >= 16];
#' fit1 <- fitdist_fr_wei(x=usecounts, param_fr_wei=param_fr_wei);
#' plot_fr_wei(fit1);
#' plot_fr_wei(fit1, ylim=c(0,0.005));
#' plot_fr_wei(fit1, scale_y=TRUE) + ggtitle("densities scaled to max height")
#' plot_fr_wei(fit1, scale_y=TRUE, highlight_range=c(100, 2200))
#'
#' @export
plot_fr_wei <- function
(fit,
 addx=NULL,
 ylim=NULL,
 xlim=NULL,
 scale_y=FALSE,
 log="x",
 type=c("gg", "base"),
 highlight_range=NULL,
 do_facet=TRUE,
 ...)
{
   ## Simple function to plot the Frechet-Weibull fit
   type <- match.arg(type);
   coef_fr_wei <- coef(fit);
   x <- unique(sort(c(addx, fit$data)));
   y1 <- dinvweibull(x=x,
      scale=coef_fr_wei["fr_scale"],
      shape=coef_fr_wei["fr_shape"],
      log=FALSE);
   y2 <- dweibull(x=x,
      scale=coef_fr_wei["wei_scale"],
      shape=coef_fr_wei["wei_shape"],
      log=FALSE);
   y12 <- y1 * coef_fr_wei["fr_weight"] +
      y2 * (1 - coef_fr_wei["fr_weight"]);
   if (scale_y) {
      y1 <- y1 / max(y1);
      y2 <- y2 / max(y2);
      y12 <- y12 / max(y12);
   }
   if (length(ylim) == 0) {
      ylim <- range(c(0, y1, y2, y12));
   }
   if (length(xlim) == 0) {
      xlim <- range(c(x, addx));
   }
   if (jamba::igrepHas("base", type)) {
      plot(x=x, y=y1, col="red", type="l",
         ylim=ylim,
         xlim=xlim,
         ...);
      lines(x=x, y=y2, col="blue", type="l");
      lines(x=x, y=y12, col="purple", type="l", lty="dotted", lwd=2);
   }
   yL <- list(y1=y1, y2=y2, y12=y12);
   dist_names <- c("Frechet", "Weibull", "Frechet-Weibull");
   dist_label <- paste0(c("Frechet", "Weibull", "Frechet-Weibull"),
      c(
         paste0("(w=",
            round(digits=2, coef_fr_wei["fr_weight"]),
            ")"),
         paste0("(w=",
            round(digits=2, 1-coef_fr_wei["fr_weight"]),
            ")"),
         ""));
   df <- data.frame(x=rep(x, 3),
      y=unlist(yL),
      dist=factor(rep(dist_names, lengths(yL)),
         levels=dist_names[c(3,1,2)]),
      dist_label=factor(rep(dist_label, lengths(yL)),
         levels=dist_label[c(3,1,2)])
   )
   if (do_facet) {
      df$type <- "Components of the Frechet-Weibull Fit";
   }

   ## Optionally highlight a range of points
   if (length(highlight_range) > 0) {
      highlight_flag <- df$x >= min(highlight_range) & df$x <= max(highlight_range);
      df$highlight <- ifelse(highlight_flag,
         "selected", "discarded");
   } else {
      df$highlight <- "discarded";
   }

   ## ggplot
   if (jamba::igrepHas("gg", type)) {
      suppressPackageStartupMessages(require(ggplot2));
      gg <- ggplot(df,
         aes(x=x,
            y=y,
            group=dist_label,
            color=dist_label,
            linetype=dist_label)) +
         geom_line(size=1) +
         colorjam::theme_jam() +
         colorjam::scale_color_jam() +
         #xlim(xlim) +
         ylim(ylim) +
         geom_vline(xintercept=coef_fr_wei[c("fr_scale", "wei_scale")]);
      if (length(highlight_range) > 0) {
         gg <- gg +
            ggforce::geom_mark_rect(#data=subset(df, highlight %in% "selected"),
               aes(filter=highlight %in% "selected" & dist %in% "Frechet-Weibull",
                  label=highlight),
               expand=unit(0.01, "mm"),
               radius=unit(0.01, "mm"));
      }
      if (do_facet) {
         gg <- gg +
            facet_wrap(~type);
      }
      if ("x" %in% log) {
         gg <- gg + scale_x_log10();
      }
      return(gg);
   }
   invisible(df);
}
