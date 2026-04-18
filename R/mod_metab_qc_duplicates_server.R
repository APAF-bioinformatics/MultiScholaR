#' @rdname mod_metab_qc_duplicates
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderUI tabsetPanel tabPanel tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
mod_metab_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        duplicate_info <- shiny::reactiveVal(NULL)
        resolution_stats <- shiny::reactiveVal(NULL)
        filter_plot <- shiny::reactiveVal(NULL)
        
        # Detect duplicates
        shiny::observeEvent(input$detect_duplicates, {
            runMetabDuplicateDetectionObserverShell(
                stateManager = workflow_data$state_manager
                , setDuplicateInfoFn = duplicate_info
            )
        })
        
        # Render duplicate summary
        output$duplicate_summary <- shiny::renderUI({
            buildMetabDuplicateSummaryUi(duplicate_info())
        })
        
        # Render duplicate tables per assay
        output$duplicate_tables <- shiny::renderUI({
            buildMetabDuplicateTablesUi(
                dupList = duplicate_info()
                , nsFn = ns
            )
        })
        
        # Render individual duplicate tables
        shiny::observe({
            registerMetabDuplicateTableRenderers(
                dupList = duplicate_info()
                , output = output
            )
        })
        
        # Resolve duplicates
        shiny::observeEvent(input$resolve_duplicates, {
            runMetabDuplicateResolutionObserver(
                workflowData = workflow_data
                , omicType = omic_type
                , output = output
                , setDuplicateInfoFn = duplicate_info
                , setResolutionStatsFn = resolution_stats
                , setFilterPlotFn = filter_plot
            )
        })
        
        # Revert
        shiny::observeEvent(input$revert_duplicates, {
            runMetabDuplicateRevertObserverShell(
                runRevertFn = function() {
                    revertMetabDuplicateResolution(
                        stateManager = workflow_data$state_manager
                    )
                }
                , output = output
                , setResolutionStatsFn = resolution_stats
                , setDuplicateInfoFn = duplicate_info
                , setFilterPlotFn = filter_plot
            )
        })

        # Render QC progress plot
        output$filter_plot <- shiny::renderPlot({
            renderMetabDuplicateFilterPlot(
                filterPlot = filter_plot
            )
        })
    })
}

