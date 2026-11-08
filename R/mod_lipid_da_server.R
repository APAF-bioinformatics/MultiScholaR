#' @rdname mod_lipid_da
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
mod_lipid_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, profile = NULL) {
    # KIT-106 W3: the shared DA frame (mod_omic_da.R). As KIT-110 D1 did for the normalisation tab, `profile` only selects
    # the omic's configuration (trx: `trxOmicProfile()`); NULL is lipidomics. The frame runs as that omic (lipid's server
    # read no `omic_type` before either).
    selected <- profile$omic_type %||% "lipidomics"
    mod_omic_da_server(id, workflow_data, experiment_paths, selected, experiment_label, config = omicGuiConfig(selected)$da)
}
