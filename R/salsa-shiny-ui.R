#' SALSA Shiny app UI
#'
#' SALSA Shiny app UI
#'
#' @family SALSA Shiny functions
#'
#' @import shiny
#' @import shinydashboard
#'
#' @export
salsaAppUI <- function
(...)
{
   # header
   header <- dashboardHeader(
      title=tagList("SALSA", icon("shoe-prints"))
   )

   # sidebar
   sidebar <- dashboardSidebar(
      sidebarMenu(
         id="tabs",
         menuItem(
            text="Guides",
            tabName="guides",
            icon=icon("info")),
         menuItem(
            text="Import Data",
            tabName="importdata",
            icon=icon("file-upload")),
         menuItem(
            text="Cell threshold",
            tabName="cellthreshold",
            icon=icon("tint")),
         menuItem(
            text="Gene threshold",
            tabName="genethreshold",
            icon=icon("dna")),
         menuItem(
            text="Export",
            tabName="export",
            icon=icon("file-download"))
      )
   );

   # Define import tab
   importTab <- fluidPage(
      shinyjs::useShinyjs(),
      fluidRow(
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Import nUMI per cell",
               status="warning",
               solidHeader=TRUE,
               width=12,
               fileInput(
                  inputId="numi_per_cell",
                  label=NULL,
                  multiple=FALSE),
               plotOutput("importCellOutput")
            )
         ),
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Import nUMI per gene",
               status="warning",
               solidHeader=TRUE,
               width=12,
               fileInput(
                  inputId="numi_per_gene",
                  label=NULL,
                  multiple=FALSE),
               plotOutput("importGeneOutput")
            )
         )
      )
   );

   # Define cellthreshold tab
   cellthresholdTab <- fluidPage(
      fluidRow(
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Calculate or Import Parameters",
               status="warning",
               solidHeader=TRUE,
               width=12,
               fileInput(
                  label="Pre-calculated nUMI per cell parameters",
                  inputId="import_cell_params",
                  multiple=FALSE),
               conditionalPanel(
                  condition="output.cell_file_uploaded == true",
                  sliderInput(
                     "cell_count_slider",
                     label="Minimum count range",
                     min=1,
                     max=2000,
                     value=c(2, 2000),
                     step=1,
                     round=TRUE,
                     pre=">="
                  ),
                  shinyBS::bsTooltip(
                     id="cell_count_slider",
                     title="Counts will be filtered for at least this many counts prior to parameter fitting.",
                     placement="top",
                     trigger="hover"
                  ),
                  actionButton("calc_cell_params",
                     label="Calculate Cell Parameters")
               )
            ),
            shinydashboard::box(
               title="Frechet summary",
               status="primary",
               solidheader=FALSE,
               width=12,
               plotlyOutput("cell_output_left",
                  height="400px")
            )
         ),
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Frechet-Weibull summary",
               status="primary",
               solidHeader=FALSE,
               width=12,
               plotlyOutput("cell_output_right",
                  height="650px")
            )
         )
      )
   );

   # Define genethreshold tab
   genethresholdTab <- fluidPage(
      fluidRow(
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Calculate or Import Parameters",
               status="warning",
               solidHeader=TRUE,
               width=12,
               fileInput(
                  label="Pre-calculated nUMI per gene parameters",
                  inputId="import_gene_params",
                  multiple=FALSE),
               conditionalPanel(
                  condition="output.gene_file_uploaded == true",
                  sliderInput(
                     "gene_count_slider",
                     label="Minimum count range",
                     min=1,
                     max=2000,
                     value=c(2, 2000),
                     step=1,
                     round=TRUE,
                     pre=">="
                  ),
                  shinyBS::bsTooltip(
                     id="gene_count_slider",
                     title="Counts will be filtered for at least this many counts prior to parameter fitting.",
                     placement="top",
                     trigger="hover"
                  ),
                  actionButton("calc_gene_params",
                     label="Calculate Gene Parameters")
               )
            ),
            shinydashboard::box(
               title="Frechet summary",
               status="primary",
               solidheader=FALSE,
               width=12,
               plotlyOutput("gene_output_left",
                  height="400px")
            )
         ),
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               title="Frechet-Weibull summary",
               status="primary",
               solidHeader=FALSE,
               width=12,
               plotlyOutput("gene_output_right",
                  height="650px")
            )
         )
      )
   );

   ## Load guidesTab via salsaAppConstants()
   salsaAppConstants();

   # dashboard body
   body <- dashboardBody(
      tabItems(
         tabItem(tabName="guides", guidesTab),
         tabItem(tabName="importdata", importTab),
         tabItem(tabName="cellthreshold", cellthresholdTab),
         tabItem(tabName="genethreshold", genethresholdTab)
         #tabItem(tabName="export", exportTab)
      )
   );

   dashboardPage(
      header,
      sidebar,
      body,
      skin="blue");
}

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
