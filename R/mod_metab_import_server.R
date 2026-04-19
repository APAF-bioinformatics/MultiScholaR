#' @rdname mod_metab_import
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput outputOptions
#' @importFrom shinyFiles shinyFileChoose parseFilePaths getVolumes
#' @importFrom logger log_info log_error log_warn
mod_metab_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  runMetabImportModuleServerPublicWrapper(
    id = id,
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = volumes
  )
}
