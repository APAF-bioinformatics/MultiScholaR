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
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::wellPanel(
            shiny::fluidRow(
                # ================================================================
                # LEFT PANEL (width = 3): All normalization controls
                # ================================================================
                if (is.null(profile)) buildLipidNormOptionsControlPanel(ns) else profile$norm_options(ns)

                # ================================================================
                # RIGHT PANEL (width = 9): QC Visualization Tabs
                # ================================================================
                # KIT-097: the profile is passed only when there is one (the lipid call is unchanged).
                , if (is.null(profile)) buildLipidNormQcTabsetPanel(ns) else buildLipidNormQcTabsetPanel(ns, profile = profile)
            )
        )
    )
}
