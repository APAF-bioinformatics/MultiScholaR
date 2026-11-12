#' Proteomics Normalization Applet Module
#'
#' A Shiny module for normalization, RUV-III batch correction, and
#' correlation-based filtering before differential-abundance analysis.
#'
#' @name normalizationAppletModule
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname normalizationAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 selectInput helpText hr numericInput checkboxInput actionButton icon br sliderInput h5 verbatimTextOutput plotOutput tabsetPanel tabPanel conditionalPanel
mod_prot_norm_ui <- function(id) {
  # KIT-120: the shared normalisation frame (mod_omic_norm.R) under proteomics' configuration (P-32 closed).
  mod_omic_norm_ui(id, omicGuiConfig("proteomics")$norm)
}
