#' Lipidomics Import Module
#'
#' A Shiny module for importing lipidomics data with vendor-format detection
#' and dynamic column mapping.
#'
#' @name mod_lipid_import
#' @param id,workflow_data,experiment_paths,volumes Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096, e.g. the transcriptomics profile); `NULL`, the default, is lipidomics.
NULL

# KIT-102 W5: a thin wrapper over the shared import frame (mod_omic_import.R) and its feature-table panel
# (mod_omic_import_feature_table.R), under `omicGuiConfig(profile$omic_type %||% "lipidomics")`. The `profile` argument
# stays for the frozen entry suite (KIT-109 removes it).

#' @rdname mod_lipid_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_lipid_import_ui <- function(id, profile = NULL) {
  mod_omic_import_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics"))
}
