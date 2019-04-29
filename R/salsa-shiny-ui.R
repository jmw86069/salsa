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
               collapsible=TRUE,
               solidHeader=TRUE,
               width=12,
               conditionalPanel(
                  condition="output.cell_count_slider == 12345678",
                  fileInput(
                     label="Pre-calculated nUMI per cell parameters",
                     inputId="import_cell_params",
                     multiple=FALSE)
               ),
               conditionalPanel(
                  condition="output.cell_file_uploaded == true",
                  shinyWidgets::prettySwitch(
                     inputId="log_cell_x",
                     label="Log-scale x-axis?",
                     fill=TRUE,
                     value=FALSE,
                     inline=TRUE
                  ),
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
                  selectInput(
                     inputId="cell_fr_wei_method",
                     label="Frechet-Weibull fit method",
                     choices=c("mle", "mge"),
                     selected="mge",
                     selectize=FALSE
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

