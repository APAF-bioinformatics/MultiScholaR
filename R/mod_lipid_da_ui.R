#' @rdname mod_lipid_da
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput textAreaInput tabsetPanel tabPanel conditionalPanel helpText
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-097: `trxOmicProfile()`): the model shown in place of the formula input, its labels and
#'   test ids, no session load and no heatmap; `NULL`, the default, is lipidomics, unchanged.
mod_lipid_da_ui <- function(id, profile = NULL) {
    # KIT-106 W3: the shared DA frame (mod_omic_da.R). As KIT-110 D1 did for the normalisation tab, `profile` only selects
    # the omic's configuration (trx: `trxOmicProfile()`); NULL is lipidomics.
    mod_omic_da_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics")$da)
}
