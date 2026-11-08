#' @title Internal Standard QC Visualization Module
#' @description A Shiny module for visualizing internal standard quality metrics.
#'              Displays IS detection, CV distributions, and intensity trends.
#'              KIT-104 W6: the exported wrappers over the shared internal-standard sub-tab (mod_omic_qc_tabs.R)
#'              under the metab configuration.
#' @name mod_metab_qc_itsd
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_metab_qc_itsd
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr textInput actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_itsd_ui <- function(id) {
    mod_omic_qc_itsd_ui(id, omicGuiConfig("metabolomics"))
}
