#' @rdname mod_lipid_da
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
mod_lipid_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    # KIT-106 W3: the shared DA frame (mod_omic_da.R) under the lipid configuration, run as lipidomics whatever
    # `omic_type` says (lipid's server read no `omic_type` before either). KIT-109: trx calls the frame itself.
    mod_omic_da_server(id, workflow_data, experiment_paths, "lipidomics", experiment_label,
        config = omicGuiConfig("lipidomics")$da)
}
