# Extracted from test-lipid-01d-qc-intensity-module-seams.R:291

# prequel ----------------------------------------------------------------------
library(testthat)
suppressPackageStartupMessages(library(shiny))
source(test_path("..", "..", "R", "mod_lipid_qc_intensity_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_qc_intensity.R"), local = environment())
ensureLipidIntensityWorkflowTestClass <- function() {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                lipid_id_column = "character"
            )
        )
    }
}

# test -------------------------------------------------------------------------
ensureLipidIntensityWorkflowTestClass()
input <- list(
        apply_filter = 7L,
        intensity_cutoff_percentile = 10L,
        proportion_below_cutoff = 0.5
    )
current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L1", "L2"),
                intensity = c(12, 9, 7),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = c("F1", "F2"),
                feature = c("F1", "F2"),
                intensity = c(4, 2),
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id",
        design_matrix = data.frame(Sample_ID = "intensity", group = "A"),
        sample_id = "Sample_ID",
        group_id = "group"
    )
filtered_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = "L1",
                intensity = 12,
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = "F1",
                feature = "F1",
                intensity = 4,
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id",
        design_matrix = data.frame(Sample_ID = "intensity", group = "A"),
        sample_id = "Sample_ID",
        group_id = "group"
    )
filter_stats_value <- NULL
filter_plot_value <- "stale-plot"
info_messages <- character()
warn_messages <- character()
notifications <- list()
removed_notifications <- character()
save_state_args <- NULL
filtering_args <- NULL
qc_plot_args <- NULL
output_env <- new.env(parent = emptyenv())
qc_plot <- structure(list(id = "qc-grid"), class = "grob")
workflow <- omicsTestOwner(current_s4, "lipid_raw_data_s4", config = list(cache = "cfg"))
wrapper <- omicsMockWrapper("runLipidReferenceIntensityFilter", list(object = filtered_s4))
result <- registerLipidIntensityApplyFilterObserver(
        input = input,
        workflowData = workflow,
        output = output_env,
        filterStatsVal = function(value) {
            filter_stats_value <<- value
        },
        filterPlotVal = function(value) {
            filter_plot_value <<- value
        },
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            expect_identical(event, 7L)
            force(handler)
            invisible(NULL)
        },
        reqFn = function(value) {
            if (is.null(value)) {
                stop("required value missing")
            }
            invisible(value)
        },
        showNotificationFn = function(message, type = NULL, id = NULL, duration = NULL) {
            notifications <<- c(notifications, list(list(
                message = message,
                type = type,
                id = id,
                duration = duration
            )))
        },
        removeNotificationFn = function(id) {
            removed_notifications <<- c(removed_notifications, id)
        },
        renderTextFn = function(expr) force(expr),
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarnFn = function(message) {
            warn_messages <<- c(warn_messages, message)
        },
        logErrorFn = function(...) {
            stop("logErrorFn should not be called on success")
        },
        updateLipidFilteringFn = function(theObject, step_name, omics_type, return_grid, overwrite) {
            qc_plot_args <<- list(
                theObject = theObject,
                step_name = step_name,
                omics_type = omics_type,
                return_grid = return_grid,
                overwrite = overwrite
            )
            qc_plot
        }
    )
expect_identical(result, input)
expect_length(wrapper$calls, 1L)
expect_identical(wrapper$calls[[1L]][[1L]], workflow$state_manager)
expect_identical(wrapper$calls[[1L]][[2L]], 10L)
expect_identical(wrapper$calls[[1L]][[3L]], 0.5)
expect_identical(wrapper$calls[[1L]][[4L]], list(cache = "cfg"))
expect_identical(qc_plot_args$theObject, filtered_s4)
expect_identical(qc_plot_args$step_name, "2_Intensity_Filtered")
expect_identical(qc_plot_args$omics_type, "lipidomics")
expect_identical(qc_plot_args$return_grid, TRUE)
expect_identical(qc_plot_args$overwrite, TRUE)
expect_equal(
        filter_stats_value,
        data.frame(
            Assay = c("AssayA", "AssayB"),
            Original = c(2, 2),
            Filtered = c(1, 1),
            Removed = c(1, 1),
            Percent_Retained = c(50, 50),
            stringsAsFactors = FALSE
        )
    )
expect_identical(filter_plot_value, qc_plot)
expect_match(output_env$filter_results, "Lipid Intensity Filter Applied Successfully", fixed = TRUE)
expect_match(output_env$filter_results, "AssayA: 2 -> 1", fixed = TRUE)
expect_match(output_env$filter_results, "AssayB: 2 -> 1", fixed = TRUE)
expect_match(output_env$filter_results, "Total: 4 -> 2 lipids (removed 2)", fixed = TRUE)
