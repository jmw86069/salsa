#' SALSA Shiny app server
#'
#' SALSA Shiny app sserver
#'
#' @param input provided by shiny
#' @param output provided by shiny
#'
#' @family SALSA Shiny functions
#'
#' @import jamba
#' @import dplyr `%>%`
#'
#' @export
salsaAppServer <- function
(input,
 output,
 session)
{

   cellthreshold_guide <- shiny::fluidPage(
      htmltools::h1("Cell Threshold",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         htmltools::tags$ul(
            htmltools::tags$li(htmltools::strong("Descriptive steps", style="color:dimgrey"),
               " go here.")
         )
      )
   );
   output$import_guide <- renderUI(import_guide);
   output$cellthreshold_guide <- renderUI(cellthreshold_guide);
   output$genethreshold_guide <- renderUI(genethreshold_guide);
   output$export_guide <- renderUI(export_guide);

   #############################################
   ## Handle cell threshold
   numi_per_cell <- reactive({
      # Require numi_per_cell file is uploaded already
      req(input$numi_per_cell);
      numi_per_cell_df <- read_numi(file=input$numi_per_cell$datapath);
      if (nrow(numi_per_cell_df) <= 1) {
         return();
      }
      if (any(c("character","factor") %in% class(numi_per_cell_df[,1])) &&
            any(c("numeric","integer") %in% class(numi_per_cell_df[,2]))) {
         shiny::showNotification("Input numi_per_cell is properly formatted",
            type="message");
      } else {
         shiny::showNotification("Input numi_per_cell must contain two columns,
            the first with 'barcodes' and the second with 'counts'.",
            duration=NULL,
            type="error");
         print(head(numi_per_cell_df));
         return();
      }
      ## update the numeric ranges allowed
      count_vector_l <- get_salsa_steps(numi_per_cell_df[,2],
         include_vector=TRUE);
      default_range <- jamba::noiseFloor(range(count_vector_l$count_vector),
         minimum=4);
      updateSliderInput(session,
         "cell_count_slider",
         min=min(count_vector_l$count_vector),
         max=max(count_vector_l$count_vector),
         value=default_range);

      return(numi_per_cell_df);
   });
   output$importCellOutput <- renderPlot({
      # Require non-empty numi_per_cell data.frame
      req(numi_per_cell());
      numi_per_cell_df <- numi_per_cell();
      gg_rank_count_cell <- rank_count_plot(numi_per_cell_df[,2]);
      gg_rank_count_cell;
   });
   output$cell_file_uploaded <- reactive({
      return(length(numi_per_cell()) > 0)
   })
   outputOptions(output,
      "cell_file_uploaded",
      suspendWhenHidden=FALSE);

   # button "calc_cell_params"
   get_cell_salsa_table <- reactive({
      ## only the button causes calculations
      input$calc_cell_params;

      ## Isolate other inputs so they do not cause immediate reaction
      cell_count_slider <- isolate(input$cell_count_slider);
      numi_per_cell_df <- isolate(numi_per_cell());
      cell_fr_wei_method <- isolate(input$cell_fr_wei_method);
      usecounts <- numi_per_cell_df[,2];

      #
      count_vector_l <- get_salsa_steps(usecounts,
         include_vector=TRUE);
      count_vector <- count_vector_l$count_vector;
      count_vector <- count_vector[count_vector >= min(cell_count_slider)
         & count_vector <= max(cell_count_slider)];
      # Currently restrict lower bound to 2
      count_vector <- count_vector[count_vector >= 2];

      ## Optionally print verbose output
      verbose <- length(getOption("verbose")) > 0 && getOption("verbose");
      if (verbose) {
         printDebug("get_cell_salsa_table(): ",
            "Starting do_salsa_steps()");
         printDebug("get_cell_salsa_table(): ",
            "length(usecounts):", length(usecounts));
         printDebug("get_cell_salsa_table(): ",
            "range(count_vector):", range(count_vector));
      }
      ## Wrap the workflow in a progress bar
      withProgress(
         message="Calculating fit parameters",
         value=0,
         {
            salsa_steps <- do_salsa_steps(x=usecounts,
               count_vector=count_vector,
               do_shiny_progress=TRUE,
               fr_wei_method=cell_fr_wei_method,
               verbose=verbose);
            salsa_table <- get_salsa_table(salsa_steps);
         }
      );
      return(salsa_table);
   });

   output$cell_output_left <- renderPlotly({
      salsa_table <- get_cell_salsa_table();
      ## Handle event data
      eventdata_r <- event_data("plotly_click", source="cell_right");
      eventdata_l <- event_data("plotly_click", source="cell_left");
      if ("data.frame" %in% class(eventdata_r)) {
         print(head(eventdata_r));
      }
      ## Plot the results
      if (input$log_cell_x) {
         xaxis_type <- "log";
      } else {
         xaxis_type <- "linear";
      }
      fr_colnames <- c("shape", "scale");
      fr_plot_list <- lapply(nameVector(fr_colnames), function(i){
         plot_ly(salsa_table,
            x=~count,
            y=as.formula(paste0("~",i)),
            name=i,
            type="scatter",
            source="cell_left",
            mode="markers+lines") %>%
         layout(
            showlegend=FALSE,
            yaxis=list(title=i),
            spikedistance=10,
            hovermode="x",
            xaxis=list(
               spikemode="across",
               spikedash="solid",
               spikecolor="grey",
               spikesnap="cursor",
               spikethickness=2,
               showspikes=TRUE,
               spikesides=TRUE,
               spikes="across",
               type=xaxis_type
            )
         ) %>%
         add_annotations(
            text=~i,
            x=0.5,
            y=1,
            yref="paper",
            xref="paper",
            yanchor="top",
            xanchor="middle",
            showarrow=FALSE
         ) %>%
         config(displaylogo=FALSE)
      });
      sp_fr <- subplot(fr_plot_list,
         nrows=length(fr_plot_list),
         shareX=TRUE);
      sp_fr;
   });
   output$cell_output_right <- renderPlotly({
      salsa_table <- get_cell_salsa_table();
      ## Plot the results
      if (input$log_cell_x) {
         xaxis_type <- "log";
      } else {
         xaxis_type <- "linear";
      }
      fr_wei_colnames <- c("fr_shape", "fr_scale",
         "fr_weight",
         "wei_shape", "wei_scale");
      fr_wei_plot_list <- lapply(nameVector(fr_wei_colnames), function(i){
         legendgroup <- gsub("^.+_", "",
            gsub("_(shape|scale)", "", i));
         plot_ly(salsa_table,
            x=~count,
            y=as.formula(paste0("~",i)),
            name=i,
            legendgroup=legendgroup,
            type="scatter",
            source="cell_right",
            mode="markers+lines") %>%
         layout(
            showlegend=FALSE,
            spikedistance=10,
            hovermode="x",
            yaxis=list(title=i),
            xaxis=list(
               spikemode="across",
               spikedash="solid",
               spikecolor="grey",
               spikesnap="cursor",
               spikethickness=2,
               showspikes=TRUE,
               spikesides=TRUE,
               spikes="across",
               type=xaxis_type
            )
         ) %>%
         add_annotations(
            text=~i,
            x=0.5,
            y=1,
            yref="paper",
            xref="paper",
            yanchor="top",
            xanchor="middle",
            showarrow=FALSE
         )
      });
      sp_fr <- subplot(fr_wei_plot_list,
         nrows=length(fr_wei_plot_list),
         shareX=TRUE);
      sp_fr;
   });

   #############################################
   ## Handle the same behavior for nUMI per gene
   numi_per_gene <- reactive({
      # Require numi_per_gene file is uploaded already
      req(input$numi_per_gene);
      numi_per_gene_df <- read_numi(input$numi_per_gene$datapath);
      if (nrow(numi_per_gene_df) <= 1) {
         return();
      }
      if (any(c("character","factor") %in% class(numi_per_gene_df[,1])) &&
            any(c("numeric","integer") %in% class(numi_per_gene_df[,2]))) {
         shiny::showNotification("Input numi_per_gene is properly formatted",
            type="message");
      } else {
         shiny::showNotification("Input numi_per_gene must contain two columns,
            the first with 'gene' and the second with 'counts'.",
            duration=NULL,
            type="error");
         return();
      }
      ## update the numeric ranges allowed
      count_vector_l <- get_salsa_steps(numi_per_gene_df[,2],
         include_vector=TRUE);
      updateSliderInput(session,
         "cell_count_slider",
         min=min(count_vector_l$count_vector),
         max=max(count_vector_l$count_vector),
         value=range(count_vector_l$count_vector));
      return(numi_per_gene_df);
   });
   output$importGeneOutput <- renderPlot({
      # Require non-empty numi_per_gene data.frame
      req(numi_per_gene());
      numi_per_gene_df <- numi_per_gene();
      gg_rank_count_gene <- rank_count_plot(numi_per_gene_df[,2]);
      gg_rank_count_gene;
   });
   output$gene_file_uploaded <- reactive({
      return(length(numi_per_gene()) > 0)
   })
   outputOptions(output,
      "gene_file_uploaded",
      suspendWhenHidden=FALSE);
   # button "calc_gene_params"
   get_gene_salsa_table <- reactive({
      ## only the button causes calculations
      input$calc_gene_params;

      ## Isolate other inputs so they do not cause immediate reaction
      gene_count_slider <- isolate(input$gene_count_slider);
      numi_per_gene_df <- isolate(numi_per_gene());
      usecounts <- numi_per_gene_df[,2];

      #
      count_vector_l <- get_salsa_steps(usecounts,
         include_vector=TRUE);
      count_vector <- count_vector_l$count_vector;
      count_vector <- count_vector[count_vector >= min(gene_count_slider)
         & count_vector <= max(gene_count_slider)];

      ## Wrap the workflow in a progress bar
      withProgress(
         message="Calculating gene fit parameters",
         value=0,
         {
            salsa_steps <- do_salsa_steps(x=usecounts,
               count_vector=count_vector,
               do_shiny_progress=TRUE,
               verbose=FALSE);
            salsa_table <- get_salsa_table(salsa_steps);
         }
      );
      return(salsa_table);
   });
   output$gene_output_left <- renderPlotly({
      salsa_table <- get_gene_salsa_table();
      ## Plot the results
      fr_colnames <- c("shape", "scale");
      fr_plot_list <- lapply(nameVector(fr_colnames), function(i){
         plot_ly(salsa_table,
            x=~count,
            y=as.formula(paste0("~",i)),
            name=i,
            type="scatter",
            mode="markers+lines") %>%
            layout(
               yaxis=list(title=i),
               spikedistance=10,
               hovermode="x",
               xaxis=list(
                  spikemode="across",
                  spikedash="solid",
                  spikecolor="grey",
                  spikesnap="cursor",
                  spikethickness=2,
                  showspikes=TRUE,
                  spikesides=TRUE,
                  spikes="across"
               )
            ) %>%
            config(displaylogo=FALSE)
      });
      sp_fr <- subplot(fr_plot_list,
         nrows=length(fr_plot_list),
         shareX=TRUE);
      sp_fr;
   });
   output$gene_output_right <- renderPlotly({
      salsa_table <- get_gene_salsa_table();
      ## Plot the results
      fr_wei_colnames <- c("fr_shape", "fr_scale",
         "fr_weight",
         "wei_shape", "wei_scale");
      fr_wei_plot_list <- lapply(nameVector(fr_wei_colnames), function(i){
         legendgroup <- gsub("^.+_", "",
            gsub("_(shape|scale)", "", i));
         plot_ly(salsa_table,
            x=~count,
            y=as.formula(paste0("~",i)),
            name=i,
            legendgroup=legendgroup,
            type="scatter",
            mode="markers+lines") %>%
            layout(
               spikedistance=10,
               hovermode="x",
               yaxis=list(title=i),
               xaxis=list(
                  spikemode="across",
                  spikedash="solid",
                  spikecolor="grey",
                  spikesnap="cursor",
                  spikethickness=2,
                  showspikes=TRUE,
                  spikesides=TRUE,
                  spikes="across"
               )
            )
      });
      sp_fr <- subplot(fr_wei_plot_list,
         nrows=length(fr_wei_plot_list),
         shareX=TRUE);
      sp_fr;
   });

}
