#' @rdname mod_metab_qc_duplicates
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderUI tabsetPanel tabPanel tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
mod_metab_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label) {
    mod_omic_qc_duplicates_server(id, workflow_data, omic_type, experiment_label, omicGuiConfig("metabolomics"))
}
