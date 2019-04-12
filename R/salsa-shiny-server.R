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
   if (!file.exists("/Users/wardjm/nchar_import_guide.txt")) {
      write.table(file="/Users/wardjm/nchar_import_guide.txt",
         data.frame(nchar_import_guide=nchar(import_guide)),
         sep="\t");
   }
   jamba::printDebug("nchar(import_guide):", nchar(import_guide));

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
      updateSliderInput(session,
         "cell_count_slider",
         min=min(count_vector_l$count_vector),
         max=max(count_vector_l$count_vector),
         value=range(count_vector_l$count_vector));

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
      usecounts <- numi_per_cell_df[,2];

      #
      count_vector_l <- get_salsa_steps(usecounts,
         include_vector=TRUE);
      count_vector <- count_vector_l$count_vector;
      count_vector <- count_vector[count_vector >= min(cell_count_slider)
         & count_vector <= max(cell_count_slider)];

      ## Wrap the workflow in a progress bar
      withProgress(
         message="Calculating fit parameters",
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
   output$cell_output_left <- renderPlotly({
      salsa_table <- get_cell_salsa_table();
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
   output$cell_output_right <- renderPlotly({
      salsa_table <- get_cell_salsa_table();
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
