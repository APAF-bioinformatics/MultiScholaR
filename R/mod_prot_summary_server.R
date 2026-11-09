#' @rdname sessionSummaryModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req renderText showNotification downloadHandler withProgress incProgress
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom purrr detect
#' @importFrom readr read_tsv
mod_prot_summary_server <- function(id, project_dirs, omic_type = "proteomics", experiment_label = NULL, workflow_data = NULL) {
  # KIT-108 W4 (D2): the shared tab (mod_omic_summary.R) under proteomics' `summary` group; its route pieces are in
  # utils_workflow_dia_summary.R.
  mod_omic_summary_server(id, project_dirs, omic_type, experiment_label, workflow_data,
    omicGuiConfig("proteomics")$summary)
}
