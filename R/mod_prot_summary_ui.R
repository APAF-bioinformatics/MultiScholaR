#' Session Summary Module
#'
#' A Shiny module for session summaries and report generation.
#'
#' @name sessionSummaryModule
#' @param id,project_dirs,omic_type,experiment_label,workflow_data Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname sessionSummaryModule
#' @export
#' @importFrom shiny NS tagList fluidPage h3 fluidRow column wellPanel h4 textInput textAreaInput actionButton icon br hr downloadButton verbatimTextOutput
mod_prot_summary_ui <- function(id) {
  # KIT-108 W4 (D2): the shared tab (mod_omic_summary.R) under proteomics' `summary` group.
  mod_omic_summary_ui(id, omicGuiConfig("proteomics")$summary)
}
