# Extracted from test-lipid-01-qc-filtering-helpers.R:740

# prequel ----------------------------------------------------------------------
library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(ggplot2))
hasLipidQcBinding <- function(name) {
    exists(name, mode = "function", inherits = TRUE)
}
skipIfMissingLipidQcBindings <- function(names) {
    missing <- names[!vapply(names, hasLipidQcBinding, logical(1))]
    skip_if(
        length(missing) > 0,
        paste("Lipid QC helper binding(s) unavailable in this ref:", paste(missing, collapse = ", "))
    )
}
hasLipidQcNamespaceBinding <- function(name) {
    exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}
skipIfMissingLipidQcNamespaceBindings <- function(names) {
    missing <- names[!vapply(names, hasLipidQcNamespaceBinding, logical(1))]
    skip_if(
        length(missing) > 0,
        paste("Lipid QC namespace binding(s) unavailable in this ref:", paste(missing, collapse = ", "))
    )
}
lipidQcNamespaceBinding <- function(name) {
    get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}
skipIfMissingFilteringProgressLipidomicsClass <- function() {
    skip_if_not(
        methods::isClass("FilteringProgressLipidomics"),
        "FilteringProgressLipidomics class unavailable in this ref"
    )
}
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
skipIfMissingLipidQcNamespaceBindings(c(
        "registerLipidIntensityApplyFilterObserver",
        "registerLipidIntensityRevertObserver"
    ))
ensureLipidIntensityWorkflowTestClass()
current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L1", "L2"),
                S1 = c(12, 9, 7),
                S2 = c(15, 10, 8),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = c("F1", "F2"),
                S1 = c(4, 2),
                S2 = c(5, 3),
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id",
        design_matrix = data.frame(
            Sample_ID = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        )
    )
filtered_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = "L1",
                S1 = 12,
                S2 = 15,
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                lipid_id = "F1",
                S1 = 4,
                S2 = 5,
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id",
        design_matrix = data.frame(
            Sample_ID = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        )
    )
input <- list(
        apply_filter = 7L,
        revert_filter = 8L,
        intensity_cutoff_percentile = 10L,
        proportion_below_cutoff = 0.5
    )
output_env <- new.env(parent = emptyenv())
filter_stats_value <- NULL
filter_plot_value <- "stale-plot"
save_state_args <- NULL
filtering_args <- NULL
qc_plot_args <- NULL
info_messages <- character()
notifications <- list()
removed_notifications <- character()
reverted_state <- NULL
qc_plot <- structure(list(id = "qc-grid"), class = "grob")
state_manager <- list(
        getState = function() current_s4,
        saveState = function(state_name, s4_data_object, config_object, description) {
            save_state_args <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
            invisible(NULL)
        },
        getHistory = function() c("initial", "lipid_intensity_filtered", "lipid_norm"),
        revertToState = function(state_name) {
            reverted_state <<- state_name
            invisible(state_name)
        }
    )
owner_workflow <- omicsTestOwner(current_s4, "lipid_raw_data_s4", config = list(cache = "cfg"))
wrapper <- omicsMockWrapper("runLipidReferenceIntensityFilter", list(object = filtered_s4))
lipidQcNamespaceBinding("registerLipidIntensityApplyFilterObserver")(
        input = input,
        workflowData = owner_workflow,
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
        logWarnFn = function(...) {
            stop("logWarnFn should not be called on apply success")
        },
        logErrorFn = function(...) {
            stop("logErrorFn should not be called on apply success")
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
expect_identical(wrapper$calls[[1L]][[1L]], owner_workflow$state_manager)
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
