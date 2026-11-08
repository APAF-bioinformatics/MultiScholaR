#' @rdname mod_lipid_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags renderImage updateSelectInput observe
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
#' @importFrom purrr imap map walk set_names
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_lipid_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL, profile = NULL) {
    # KIT-105 W4: the shared normalisation frame (mod_omic_norm.R). KIT-110 D1: `profile` only selects the omic's
    # configuration (trx: `trxOmicProfile()`); NULL is lipidomics.
    mod_omic_norm_server(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab,
        config = omicGuiConfig(profile$omic_type %||% "lipidomics")$norm)
}
