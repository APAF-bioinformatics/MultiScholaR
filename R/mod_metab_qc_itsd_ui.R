#' @rdname mod_metab_qc_itsd
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr textInput actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_itsd_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Internal Standards"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Internal Standard QC")
                    , shiny::p("Visualize internal standard (IS) performance across samples. Good IS should have low CV and consistent intensities.")
                    , shiny::hr()
                    
                    # IS pattern input
                    , shiny::textInput(
                        ns("is_pattern")
                        , "Internal Standard Regex Pattern"
                        , value = ""
                        , placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"
                    )
                    , shiny::helpText("Regular expression to identify internal standards in metabolite IDs. Leave empty to use pattern from S4 object.")
                    
                    , shiny::hr()
                    , shiny::actionButton(
                        ns("analyze_is")
                        , "Analyze Internal Standards"
                        , class = "btn-primary"
                        , width = "100%"
                    )
                    
                    , shiny::hr()
                    , shiny::h5("IS Detection Summary")
                    , shiny::uiOutput(ns("is_summary"))
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("is_results"))
                , shiny::br()
                # Per-assay IS visualization tabs
                , shiny::uiOutput(ns("is_viz_tabs"))
            )
        )
    )
}

