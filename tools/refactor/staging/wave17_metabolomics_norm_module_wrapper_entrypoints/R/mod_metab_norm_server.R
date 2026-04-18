#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags renderImage updateSelectInput observe
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
#' @importFrom purrr imap map walk set_names
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_metab_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        logger::log_info("=== METABOLOMICS NORMALIZATION MODULE STARTED ===")
        logger::log_info(paste("Module ID:", id))

        # ================================================================
        # Local Reactive Values
        # ================================================================
        norm_data <- shiny::reactiveValues(
            # --- Assay information ---
            assay_names = NULL
            , itsd_selections = list()  # Per-assay ITSD selections

            # --- Processing stages ---
            , pre_norm_qc_generated = FALSE
            , normalization_complete = FALSE
            , ruv_complete = FALSE
            , correlation_filtering_complete = FALSE

            # --- QC plot file paths (keyed by assay_name and stage) ---
            , qc_plot_paths = list()

            # --- S4 Objects at different stages ---
            , post_filter_obj = NULL
            , post_itsd_obj = NULL
            , post_log2_obj = NULL
            , post_norm_obj = NULL
            , ruv_corrected_obj = NULL
            , correlation_filtered_obj = NULL

            # --- RUV optimization results (per-assay) ---
            , ruv_optimization_results = list()

            # --- Correlation filtering ---
            , correlation_results = list()

            # --- Plot refresh trigger ---
            , plot_refresh_trigger = 0

            # --- Normalization log ---
            , normalization_log = character(0)
        )

        # ================================================================
        # Helper Functions
        # ================================================================

        # Add message to normalization log
        add_log <- function(message) {
            appendMetabNormNormalizationLog(
                normData = norm_data
                , message = message
            )
        }

        # ================================================================
        # Initialize Assay Names from S4 Object
        # ================================================================
        shiny::observe({
            initializeMetabNormAssayNames(
                stateManager = workflow_data$state_manager
                , normData = norm_data
            )
        })

        # ================================================================
        # Auto-trigger Pre-Normalization QC when tab is selected
        # ================================================================
        if (!is.null(selected_tab)) {
            shiny::observeEvent(selected_tab(), {
                runMetabNormAutoPreNormalizationQcObserverShell(
                    selectedTab = selected_tab()
                    , workflowData = workflow_data
                    , experimentPaths = experiment_paths
                    , normData = norm_data
                    , colorVariable = input$color_variable
                    , shapeVariable = input$shape_variable
                    , addLogFn = add_log
                )
            }, ignoreInit = FALSE)
        }

        # ================================================================
        # Update Plot Aesthetic Choices from Design Matrix
        # ================================================================
        shiny::observe({
            updateMetabNormDesignDrivenChoices(
                session = session
                , designMatrix = workflow_data$design_matrix
            )
        })

        # ================================================================
        # Render Normalization Log
        # ================================================================
        output$norm_log <- renderMetabNormNormalizationLog(
            normData = norm_data
        )

        # ================================================================
        # Dynamic UI: ITSD Selection Tables (per-assay)
        # ================================================================
        output$itsd_selection_ui <- renderMetabNormItsdSelectionUi(
            normData = norm_data
            , ns = ns
        )

        # ================================================================
        # Dynamic UI: RUV QC Plots (per-assay stacked)
        # ================================================================
        output$ruv_qc_ui <- renderMetabNormRuvQcUi(
            normData = norm_data,
            ns = ns
        )

        # ================================================================
        # STATIC Output Bindings: QC Plot Images (24 total)
        # Bound at server startup - no race condition
        # ================================================================
        runMetabNormQcImageBindingShell(
            output = output,
            normData = norm_data,
            qcDir = experiment_paths$metabolite_qc_dir
        )
        
        # ================================================================
        # STATIC Output Bindings: Assay Labels (8 total)
        # ================================================================
        runMetabNormAssayLabelBindingShell(
            output = output,
            getAssayNamesFn = function() norm_data$assay_names
        )

        # ================================================================
        # Render ITSD Selection Tables (per-assay)
        # ================================================================
        shiny::observe({
            runMetabNormItsdSelectionTableObserverShell(
                normData = norm_data,
                workflowData = workflow_data,
                output = output
            )
        })

        # ================================================================
        # Track ITSD Selections from DT
        # ================================================================
        shiny::observe({
            runMetabNormItsdSelectionTrackingObserverShell(
                normData = norm_data,
                input = input
            )
        })

        # ================================================================
        # Main Normalization Pipeline
        # ================================================================
        shiny::observeEvent(input$run_normalization, {
            runMetabNormNormalizationObserverWrapper(
                workflowData = workflow_data
                , input = input
                , experimentPaths = experiment_paths
                , omicType = omic_type
                , normData = norm_data
                , addLogFn = add_log
                , showNotificationFn = shiny::showNotification
                , reqFn = shiny::req
                , withProgressFn = shiny::withProgress
            )
        })

        # ================================================================
        # Reset to Pre-Normalization
        # ================================================================
        shiny::observeEvent(input$reset_normalization, {
            runMetabNormResetNormalizationObserverWrapper(
                workflowData = workflow_data
                , normData = norm_data
                , addLogFn = add_log
                , showNotificationFn = shiny::showNotification
                , reqFn = shiny::req
            )
        })

        # ================================================================
        # Render RUV Cancor Plots (per-assay)
        # ================================================================
        shiny::observe({
            runMetabNormRuvBindingObserverShell(
                normData = norm_data,
                output = output
            )
        })

        # ================================================================
        # Correlation Filtering
        # ================================================================
        shiny::observeEvent(input$apply_correlation_filter, {
            runMetabNormApplyCorrelationObserverWrapper(
                workflowData = workflow_data
                , input = input
                , normData = norm_data
                , addLogFn = add_log
                , showNotificationFn = shiny::showNotification
                , removeNotificationFn = shiny::removeNotification
                , reqFn = shiny::req
            )
        })

        # ================================================================
        # Skip Correlation Filtering
        # ================================================================
        shiny::observeEvent(input$skip_correlation_filter, {
            runMetabNormSkipCorrelationObserverWrapper(
                workflowData = workflow_data
                , normData = norm_data
                , addLogFn = add_log
                , showNotificationFn = shiny::showNotification
                , reqFn = shiny::req
            )
        })

        # ================================================================
        # Correlation Filter Summary
        # ================================================================
        output$correlation_filter_summary <- renderMetabNormCorrelationFilterSummary(
            normData = norm_data
        )

        # ================================================================
        # Final QC Plot
        # ================================================================
        output$final_qc_plot <- renderMetabNormFinalQcPlot(
            normData = norm_data
            , colorVariableFn = function() input$color_variable
            , shapeVariableFn = function() input$shape_variable
        )

        # ================================================================
        # Export Session (Comprehensive - matches proteomics pattern)
        # ================================================================
        shiny::observeEvent(input$export_session, {
            runMetabNormExportSessionObserverWrapper(
                workflowData = workflow_data
                , input = input
                , normData = norm_data
                , experimentPaths = experiment_paths
                , experimentLabel = experiment_label
                , addLogFn = add_log
                , logInfoFn = logger::log_info
                , reqFn = shiny::req
            )
        })
    })
}

