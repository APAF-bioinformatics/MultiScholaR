#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags renderImage updateSelectInput observe
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
#' @importFrom purrr imap map walk set_names
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_metab_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    # KIT-105 W5: the shared normalisation frame (mod_omic_norm.R) under the metab configuration.
    mod_omic_norm_server(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab,
        config = omicGuiConfig("metabolomics")$norm)
}
