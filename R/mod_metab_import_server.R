#' @rdname mod_metab_import
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput outputOptions
#' @importFrom shinyFiles shinyFileChoose parseFilePaths getVolumes
#' @importFrom logger log_info log_error log_warn
mod_metab_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  # KIT-102 W5: the shared import frame under the metabolomics configuration. Its completion opens the Design step.
  mod_omic_import_server(id, workflow_data, experiment_paths, volumes, config = omicGuiConfig("metabolomics"))
}
