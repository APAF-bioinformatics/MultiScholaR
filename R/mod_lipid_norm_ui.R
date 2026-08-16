#' Lipidomics Normalization Module
#'
#' A Shiny module for multi-step lipidomics normalization, harmonized with the
#' proteomics workflow.
#'
#' @name mod_lipid_norm
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_lipid_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
mod_lipid_norm_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::wellPanel(
            shiny::fluidRow(
                # ================================================================
                # LEFT PANEL (width = 3): All normalization controls
                # ================================================================
                buildLipidNormOptionsControlPanel(ns)

                # ================================================================
                # RIGHT PANEL (width = 9): QC Visualization Tabs
                # ================================================================
                , buildLipidNormQcTabsetPanel(ns)
            )
        )
    )
}
