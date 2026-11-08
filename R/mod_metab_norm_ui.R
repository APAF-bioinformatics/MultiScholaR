#' Metabolomics Normalization Module
#'
#' A Shiny module for multi-step metabolomics normalization, harmonized with
#' the proteomics workflow.
#'
#' @name mod_metab_norm
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
mod_metab_norm_ui <- function(id) {
    # KIT-105 W5: the shared normalisation frame (mod_omic_norm.R) under the metab configuration.
    mod_omic_norm_ui(id, omicGuiConfig("metabolomics")$norm)
}
