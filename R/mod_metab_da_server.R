#' @rdname mod_metab_da
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
mod_metab_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    # KIT-106 W4: the shared DA frame (mod_omic_da.R) under metabolomics' `da` configuration. The frame runs as
    # metabolomics (metab's server read no `omic_type` before either).
    mod_omic_da_server(id, workflow_data, experiment_paths, "metabolomics", experiment_label,
        config = omicGuiConfig("metabolomics")$da)
}
