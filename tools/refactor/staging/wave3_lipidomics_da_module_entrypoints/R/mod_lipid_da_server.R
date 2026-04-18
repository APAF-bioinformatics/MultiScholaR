#' @rdname mod_lipid_da
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler observe renderPrint
#' @importFrom DT renderDT datatable formatRound formatStyle styleEqual
#' @importFrom logger log_info log_error log_warn
mod_lipid_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ================================================================
        # LOCAL REACTIVE VALUES
        # ================================================================
        da_data <- initializeLipidDaServerState()

        # ================================================================
        # CONTRASTS DISPLAY
        # ================================================================
        registerLipidDaPrimaryTextOutputs(
            output = output,
            workflowData = workflow_data,
            daData = da_data
        )

        # ================================================================
        # LOAD FILTERED SESSION
        # ================================================================
        shiny::observeEvent(input$load_filtered_session, {
            bootstrapLipidDaLoadFilteredSession(
                experimentPaths = experiment_paths,
                workflowData = workflow_data,
                daData = da_data,
                session = session
            )
        })

        # ================================================================
        # RUN DA ANALYSIS
        # ================================================================
        shiny::observeEvent(input$run_da_analysis, {
            bootstrapLipidDaRunAnalysis(
                daData = da_data,
                workflowData = workflow_data,
                formulaString = input$formula_string,
                daQValThresh = input$da_q_val_thresh,
                treatLfcCutoff = input$treat_lfc_cutoff,
                session = session,
                experimentPaths = experiment_paths
            )
        })

        # ================================================================
        # VOLCANO PLOT - GLIMMA
        # ================================================================
        # ================================================================
        # WARNING BANNER
        # ================================================================
        registerLipidDaHeatmapWarningOutput(
            output = output,
            daData = da_data
        )

        registerLipidDaVolcanoGlimmaOutput(
            output = output,
            input = input,
            daData = da_data
        )

        # ================================================================
        # VOLCANO PLOT - STATIC
        # ================================================================
        registerLipidDaVolcanoStaticOutput(
            output = output,
            input = input,
            daData = da_data
        )

        # ================================================================
        # HEATMAP
        # ================================================================
        registerLipidDaHeatmapPlotOutput(
            output = output,
            input = input,
            daData = da_data
        )

        # Cluster Summary Output
        registerLipidDaClusterSummaryOutput(
            output = output,
            input = input,
            daData = da_data
        )

        # Save Heatmap Observer
        registerLipidDaSaveHeatmapObserver(
            input = input,
            daData = da_data,
            experimentPaths = experiment_paths
        )

        # ================================================================
        # DA RESULTS TABLE
        # ================================================================
        registerLipidDaSummaryStatsOutput(
            output = output,
            input = input,
            daData = da_data
        )

        registerLipidDaResultsTableOutput(
            output = output,
            input = input,
            daData = da_data
        )

        # ================================================================
        # DOWNLOAD HANDLER
        # ================================================================
        registerLipidDaResultsDownloadOutput(
            output = output,
            daData = da_data
        )
    })
}

