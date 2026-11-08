#' Lipidomics Normalization Module
#'
#' A Shiny module for multi-step lipidomics normalization, harmonized with the
#' proteomics workflow.
#'
#' @name mod_lipid_norm
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-097: `trxOmicProfile()`): its options panel, QC panels and server; `NULL`, the
#'   default, is lipidomics, unchanged.
NULL

#' @rdname mod_lipid_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
mod_lipid_norm_ui <- function(id, profile = NULL) {
    # KIT-105 W4: the shared normalisation frame (mod_omic_norm.R). KIT-110 D1: `profile` only selects the omic's
    # configuration (trx: `trxOmicProfile()`); NULL is lipidomics.
    mod_omic_norm_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics")$norm)
}
