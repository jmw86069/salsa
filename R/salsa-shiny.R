#' Launch SALSA R-Shiny app
#'
#' Launch SALSA R-Shiny app
#'
#' This function launches the SALSA workflow as an
#' R-shiny app.
#'
#' @family SALSA Shiny functions
#'
#' @param ... additional arguments are ignored.
#'
#' @examples
#' # Note: disabled for web page examples
#' # launchSalsaApp();
#'
#' @export
launchSalsaApp <- function
(...)
{
   shiny::shinyApp(ui=salsaAppUI,
      server=salsaAppServer,
      onStart=salsaAppConstants);
}
