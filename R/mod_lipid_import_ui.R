#' @rdname mod_lipid_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_lipid_import_ui <- function(id) {
  ns <- shiny::NS(id)

  # Check if shinyFiles is available
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)

  buildLipidImportUiShell(
    ns = ns,
    useShinyFiles = use_shiny_files
  )
}

