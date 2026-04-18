# Keep the auto pre-normalization QC helper top-level so a later exact-source
# wave can move it without reopening the selected-tab observer body.
generateMetabNormPreNormalizationQc <- function(
    workflowData,
    experimentPaths,
    normData,
    getPlotAestheticsFn,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    generateMetabQcPlotsFn = generateMetabQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    logInfoFn("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")

    reqFn(workflowData$state_manager)
    current_s4 <- workflowData$state_manager$getState()

    if (is.null(current_s4)) {
        logWarnFn("No S4 object available for QC plot generation")
        return()
    }

    if (inherits(current_s4, "MetaboliteAssayData")) {
        detected_assays <- names(current_s4@metabolite_data)
        if (length(detected_assays) > 0 && is.null(normData$assay_names)) {
            normData$assay_names <- detected_assays
            logInfoFn(paste("Set assay names:", paste(detected_assays, collapse = ", ")))
        }
    }

    aesthetics <- getPlotAestheticsFn()

    tryCatch({
        generateMetabQcPlotsFn(
            theObject = current_s4
            , experiment_paths = experimentPaths
            , stage = "post_filter"
            , grouping_variable = aesthetics$color_var
            , shape_variable = aesthetics$shape_var
        )

        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        normData$pre_norm_qc_generated <- TRUE
        logInfoFn("Pre-normalization QC plots generated successfully")

    }, error = function(e) {
        logErrorFn(paste("Error generating pre-normalization QC:", e$message))
        addLogFn(paste("Error generating Pre-QC:", e$message))
    })
}

