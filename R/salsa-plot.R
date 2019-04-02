
#' Rank-Count plot
#'
#' Rank-Count plot
#'
#' This function creates a classic "knee plot" showing
#' `"Rank"` on the x-axis, and `"Count"` on the y-axis,
#' typically with log-scaled x- and y-axes.
#'
#' @param x numeric or integer vector of counts
#' @param ... additional arguments are ignored
#'
#' @export
rank_count_plot <- function
(x,
 UMIname="UMIname",
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
   if (1 == 2) {
      yvals <- as.numeric(0 + CellData[,"rank"]);
      #ylab <- "UMI rank";
      xlab <- "Count";
      xvals <- as.numeric((0 + CellData[,"Count"]));

      ## data.frame with counts by rank
      rcdf <- data.frame(x=xvals,
         y=yvals,
         name=rep(UMIname, length.out=length(xvals)),
         type="Rank-Count");
   }

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
         ggforce::geom_mark_rect(data=subset(df, highlight %in% "selected"),
            aes(#filter=highlight=="selected",
               label=highlight),
            label.buffer=unit(0, "mm"),
            expand=0,
            radius=0);
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

