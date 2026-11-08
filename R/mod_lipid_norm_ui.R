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
    # KIT-105 W4: the shared normalisation frame (mod_omic_norm.R) under the lipid configuration.
    norm <- omicGuiConfig("lipidomics")$norm
    if (is.null(profile)) return(mod_omic_norm_ui(id, norm))
    # KIT-097 (KIT-105 D4, until KIT-110): the profile's options panel, and the frame's tabset with the profile's tabs, stage
    # captions and RUV texts, one assay row.
    norm$panels <- profile$norm_options
    norm$plot_tabs <- profile$norm_tabs %||% norm$plot_tabs
    norm$stage_labels <- profile$norm_stage_labels %||% norm$stage_labels
    norm$assay_slots <- 1L
    norm$ruv_heading <- profile$norm_ruv_heading %||% norm$ruv_heading
    norm$ruv_help <- profile$norm_ruv_help %||% norm$ruv_help
    mod_omic_norm_ui(id, norm)
}
