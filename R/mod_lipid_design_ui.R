#' @rdname lipidomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 h4 p conditionalPanel div icon tags HTML fluidRow column actionButton
#' @importFrom DT DTOutput renderDT
mod_lipid_design_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::fluidRow(
                        shiny::column(8
                            , shiny::h3("Design Matrix Builder")
                        )
                        , shiny::column(4
                            , shiny::actionButton(ns("show_import_modal"), "Import Existing Design"
                                , icon = shiny::icon("folder-open"), class = "btn-info pull-right")
                        )
                    )
                    , shiny::p("Use the tools below to define your experimental groups and contrasts for differential analysis. All assays (Pos/Neg modes) share the same design matrix.")

                    # Multi-assay info banner
                    , shiny::div(
                        class = "alert alert-info"
                        , style = "margin-top: 10px;"
                        , shiny::icon("info-circle")
                        , shiny::tags$strong(" Multi-Assay Workflow: ")
                        , "Design applies uniformly across all assay modes (LCMS_Pos, LCMS_Neg, GCMS, etc.). "
                        , "Sample assignments and contrasts are shared."
                    )

                    # Conditional panel: show builder when data is available
                    , shiny::conditionalPanel(
                        condition = paste0("output['", ns("data_available"), "']")
                        , if (exists("mod_lipid_design_builder_ui")) {
                            mod_lipid_design_builder_ui(ns("builder"))
                        } else {
                            shiny::div("Design builder module not loaded")
                        }

                        # Saved results preview
                        , shiny::hr()
                        , shiny::h3("Saved Results Preview")
                        , shiny::p("This section shows the design matrix and contrasts saved to the workflow.")
                        , shiny::conditionalPanel(
                            condition = paste0("output['", ns("design_matrix_exists"), "']")
                            , shiny::wellPanel(
                                shiny::h4("Current Design Matrix")
                                , DT::DTOutput(ns("design_matrix_preview"))
                                , shiny::br()
                                , shiny::h4("Defined Contrasts")
                                , DT::DTOutput(ns("contrasts_preview"))
                                , shiny::br()
                                , shiny::h4("Assays Included")
                                , shiny::verbatimTextOutput(ns("assays_preview"))
                            )
                        )
                    )

                    # Conditional panel: message if no data
                    , shiny::conditionalPanel(
                        condition = paste0("!output['", ns("data_available"), "']")
                        , shiny::div(
                            class = "alert alert-info"
                            , shiny::icon("info-circle")
                            , " Please complete the 'Import' step first. The builder will appear once data is available."
                        )
                    )
                )
            )
        )
    )
}

