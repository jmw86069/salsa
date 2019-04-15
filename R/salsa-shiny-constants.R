
#' SALSA Shiny app constants
#'
#' SALSA Shiny app constants
#'
#' This function is intended to define constant values used in
#' the creation of the SALSA shiny UI.
#'
#' This function is intended to be called only from within
#' an R-shiny app, since it defines several variables in the
#' parent (global) environment.
#'
#' @family SALSA Shiny functions
#'
#' @import shiny
#' @import shinydashboard
#' @import htmltools
#'
#' @export
salsaAppConstants <- function
(...)
{
   #
   # guides
   # define guides tab
   guidesTab <<- fluidPage(
      tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
      sidebarLayout(
         mainPanel(
            width=8,
            tabBox(
               width=12,
               tabPanel(
                  title="Import Data",
                  uiOutput("import_guide")),
               tabPanel(
                  title="Cell Threshold",
                  uiOutput("cellthreshold_guide")),
               tabPanel(
                  title="Gene Threshold",
                  verbatimTextOutput("summary")),
               tabPanel(
                  title="Export",
                  uiOutput("export_guide")))),
         sidebarPanel(
            width=4,
            "Single-cell Amalgamation by Latent Semantic Analysis (SALSA)
            is intended for analysis of single cell RNA-seq (scRNA-seq)
            data. The workflow begins by defining the useful subset of
            cells and genes to use for downstream analysis.",
            tags$ul(
               tags$li(
                  strong(style="color:firebrick",
                     "The preprint to this project is accesible on"), strong("bioRxiv"),
                  strong(style="color:firebrick", "at"), a("Lozoya et al., 2019:",
                     em("Patterns, Profiles, and Parsimony: dissecting transcriptional signatures from minimal single-cell RNA-seq output with SALSA."),
                     href="https://doi.org/10.1101/551762")
               )
            )
         )
   )
      );

   # import guide
   import_guide <<- fluidPage(
      h1("Import Data",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(strong("Import nUMI per Cell", style="color:navy"),
               "- tab-delimited file with the number of UMI per cell",
               tags$ul(
                  tags$li(
                     strong("barcodes"),
                     " in the first column"
                  ),
                  tags$li(
                     strong("counts"),
                     " in the second column"
                  ),
                  tags$li("proceed to ",
                     strong("Cell Threshold", style="color:dimgrey"),
                     " to continue the analysis.")
               )
            ),
            tags$li(strong("Import nUMI per Gene", style="color:navy"),
               "- tab-delimited file with the number of UMI per gene",
               tags$ul(
                  tags$li(
                     strong("genes"),
                     " in the first column"
                  ),
                  tags$li(
                     strong("counts"),
                     " in the second column"
                  ),
                  tags$li("proceed to ",
                     strong("Gene Threshold", style="color:dimgrey"),
                     " to continue the analysis.")
               )
            )
         )
      )
   );
   cellthreshold_guide <<- fluidPage(
      h1("Cell Threshold",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(strong("Descriptive steps", style="color:dimgrey"),
               " go here.")
         )
      )
   );
   genethreshold_guide <<- fluidPage(
      h1("Gene Threshold",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(strong("Descriptive steps", style="color:dimgrey"),
               " go here.")
         )
      )
   );
   export_guide <<- fluidPage(
      h1("Export Results",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(strong("Descriptive steps", style="color:dimgrey"),
               " go here.")
         )
      )
   );
   invisible(NULL);
}
