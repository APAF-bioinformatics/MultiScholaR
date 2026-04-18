#' @rdname mod_metab_qc_itsd
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI renderPlot tabsetPanel tabPanel tags
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_line labs theme_minimal theme element_text coord_flip scale_fill_brewer facet_wrap geom_hline
#' @importFrom logger log_info log_error log_warn
mod_metab_qc_itsd_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        runMetabQcItsdServerBody(
            input = input,
            output = output,
            session = session,
            workflowData = workflow_data,
            omicType = omic_type,
            experimentLabel = experiment_label
        )
    })
}

