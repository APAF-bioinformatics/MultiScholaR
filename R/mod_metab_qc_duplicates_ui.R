#' Metabolite Duplicate Resolution Module
#'
#' A Shiny module for identifying and resolving duplicate metabolite features.
#' KIT-104 W6: the exported wrappers over the shared duplicate sub-tab (mod_omic_qc_tabs.R) under the metab
#' configuration.
#'
#' @name mod_metab_qc_duplicates
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_metab_qc_duplicates
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_duplicates_ui <- function(id) {
    mod_omic_qc_duplicates_ui(id, omicGuiConfig("metabolomics"))
}
