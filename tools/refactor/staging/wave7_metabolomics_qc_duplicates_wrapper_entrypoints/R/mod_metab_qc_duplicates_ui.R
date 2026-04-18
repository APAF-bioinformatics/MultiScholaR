#' @rdname mod_metab_qc_duplicates
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_duplicates_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Duplicate Resolution"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Duplicate Feature Resolution")
                    , shiny::p("Identify and resolve duplicate metabolite IDs within each assay. Duplicates are resolved by keeping the feature with the highest average intensity across samples.")
                    , shiny::hr()

                    # Detection summary
                    , shiny::h5("Detection Summary")
                    , shiny::uiOutput(ns("duplicate_summary"))

                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("detect_duplicates")
                            , "Detect Duplicates"
                            , class = "btn-info"
                            , width = "100%"
                        )
                    )
                    , shiny::br()
                    , shiny::div(
                        shiny::actionButton(
                            ns("resolve_duplicates")
                            , "Resolve Duplicates"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_duplicates")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("resolution_results"))
                , shiny::br()
                # Per-assay duplicate tables
                , shiny::uiOutput(ns("duplicate_tables"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

