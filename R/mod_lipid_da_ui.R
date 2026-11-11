#' @rdname mod_lipid_da
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput textAreaInput tabsetPanel tabPanel conditionalPanel helpText
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
mod_lipid_da_ui <- function(id) {
    # KIT-106 W3: the shared DA frame (mod_omic_da.R) under the lipid configuration (KIT-109: trx calls the frame itself).
    mod_omic_da_ui(id, omicGuiConfig("lipidomics")$da)
}
