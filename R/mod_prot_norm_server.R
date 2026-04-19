#' @rdname normalizationAppletModule 
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification renderPlot renderText withProgress incProgress observe updateSelectInput
#' @importFrom readr write_tsv
#' @importFrom grid grid.draw grid.newpage pushViewport popViewport viewport grid.layout rasterGrob
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom vroom vroom_write
#' @importFrom dplyr filter select mutate distinct
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_prot_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  runProtNormModuleServerPublicWrapper(
    id = id,
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    omic_type = omic_type,
    experiment_label = experiment_label,
    selected_tab = selected_tab
  )
}
