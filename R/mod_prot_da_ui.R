#' Differential Abundance Applet Module
#'
#' A Shiny module for proteomics differential-abundance analysis and result
#' visualization across user-defined contrasts.
#'
#' @name differentialAbundanceAppletModule
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname differentialAbundanceAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 h5 textAreaInput helpText hr numericInput actionButton icon verbatimTextOutput br downloadButton
mod_prot_da_ui <- function(id) {
  # KIT-106 W5 (D5): the shared DA tab's UI from proteomics' `da` group (its layout data, mod_omic_da_panels.R).
  mod_omic_da_ui(id, omicGuiConfig("proteomics")$da)
}

#' Proteomics' archive chooser, the `archive` plugin of its DA group (mod_prot_da_ui.R:81-88 at fe97aac1): an optional
#' external native archive or historical RDS; no selection uses the latest.
#' @param ns The module's namespace function.
#' @noRd
prot_da_session_archive <- function(ns) {
  if (requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()) {
    shinyFiles::shinyFilesButton(ns("session_archive"), "Choose saved session (optional)",
      "Choose DIA session", multiple = FALSE)
  } else {
    shiny::fileInput(ns("session_archive_standard"), "Saved session (optional):",
      accept = c(".msdia.tar", ".rds", ".RDS"))
  }
}
