#' Metabolomics Import Module
#'
#' A Shiny module for importing metabolomics data with vendor-format detection
#' and dynamic column mapping.
#'
#' @name mod_metab_import
#' @param id,workflow_data,experiment_paths,volumes Runtime inputs used by this function; see the usage section for accepted values.
NULL

# KIT-102 W5: a thin wrapper over the shared import frame (mod_omic_import.R) and its feature-table panel
# (mod_omic_import_feature_table.R), under `omicGuiConfig("metabolomics")`. Each id, label and message stays metab's; the
# layout is the shared panel's.

#' @rdname mod_metab_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_metab_import_ui <- function(id) {
  mod_omic_import_ui(id, omicGuiConfig("metabolomics"))
}
