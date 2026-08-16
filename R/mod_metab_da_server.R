#' @rdname mod_metab_da
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
mod_metab_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    runMetabDaServerEntry(
        id = id,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        omicType = omic_type,
        experimentLabel = experiment_label
    )
}

