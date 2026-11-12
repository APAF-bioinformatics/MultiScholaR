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
  # KIT-120: the shared normalisation frame (mod_omic_norm.R) under proteomics' configuration; the tab's actions and
  # outputs are proteomics' own server on it (protNormFrameServer, mod_prot_norm_observer_helpers.R).
  mod_omic_norm_server(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab,
    config = omicGuiConfig("proteomics")$norm)
}
