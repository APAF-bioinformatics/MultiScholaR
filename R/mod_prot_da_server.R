#' @rdname differentialAbundanceAppletModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification renderPlot renderText withProgress incProgress
#' @importFrom purrr detect map set_names
#' @importFrom dplyr filter select mutate bind_rows slice
#' @importFrom stringr str_detect str_extract
#' @importFrom DT renderDT datatable formatRound
mod_prot_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  # KIT-106 W5 (D5): the shared DA tab with proteomics' `da` group; its init, load and run handlers are the group's
  # plugins and keep its reactive values, which the tab returns.
  mod_omic_da_server(id, workflow_data, experiment_paths, omic_type, experiment_label, omicGuiConfig("proteomics")$da,
    selected_tab = selected_tab)
}
