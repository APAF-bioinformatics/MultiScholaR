#' @title create_s4_from_protein_data_ui
#'
#' @description UI component for creating a ProteinQuantitativeData S4 object
#' from protein-level data. This provides a user-facing button to trigger
#' the process for workflows like TMT.
#'
#' @param id Module ID
#'
#' @export
#' @import shiny
create_s4_from_protein_data_ui <- function(id) {
  ns <- NS(id)
  
  shiny::fluidPage(
    shiny::wellPanel(
      shiny::h4("Finalise Protein Data"),
      shiny::p("This workflow starts with protein-level data. Click the button below to format the data into an S4 object for downstream analysis."),
      shiny::actionButton(ns("create_protein_s4"), "Create Protein S4 Object", class = "btn-primary"),
      shiny::hr(),
      shiny::h5("Results"),
      shiny::verbatimTextOutput(ns("s4_creation_results"))
    )
  )
}
