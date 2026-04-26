# fidelity-coverage-compare: shared
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

test_that("prepareLipidFilteringContext resolves slot defaults and assay naming", {
    skipIfMissingLipidQcBindings(c(
        "createLipidomicsAssayData",
        "prepareLipidFilteringContext"
    ))

    assay_data <- data.frame(
        lipid_id = c("L1", "L2"),
        S1 = c(10, 20),
        S2 = c(30, 40),
        check.names = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group",
        internal_standard_regex = "^IS_"
    )

    context <- prepareLipidFilteringContext(lipid_object)

    expect_named(
        context,
        c(
            "assay_list",
            "assay_names",
            "design_matrix",
            "group_id_col",
            "sample_id_col",
            "lipid_id_col",
            "is_pattern",
            "sample_columns"
        )
    )
    expect_identical(context$assay_names, "Assay_1")
    expect_identical(names(context$assay_list), "Assay_1")
    expect_identical(context$group_id_col, "group")
    expect_identical(context$sample_id_col, "Run")
    expect_identical(context$lipid_id_col, "lipid_id")
    expect_identical(context$is_pattern, "^IS_")
    expect_identical(context$design_matrix$Run, c("S1", "S2"))
    expect_identical(context$sample_columns, c("S1", "S2"))
})

test_that("calculateLipidFilteringAssayMetrics returns the expected assay summary contract", {
    skipIfMissingLipidQcBindings(c(
        "calculateLipidFilteringAssayMetrics",
        "countUniqueLipids",
        "countLipidsPerSample",
        "calculateLipidMissingness",
        "calculateLipidSumIntensityPerSample",
        "calculateLipidCVs",
        "getLipidInternalStandardMetrics"
    ))

    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        annotation = c("A", "B", "C"),
        S1 = c(10, 0, 5),
        S2 = c(10, NA, 15),
        S3 = c(20, 0, 0),
        check.names = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3"),
        group = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )

    metrics <- calculateLipidFilteringAssayMetrics(
        current_assay_data = assay_data,
        current_assay_name = "Assay1",
        design_matrix = design_matrix,
        group_id_col = "group",
        sample_id_col = "Run",
        lipid_id_col = "lipid_id",
        is_pattern = NA_character_,
        sample_columns = c("S1", "S2", "S3")
    )

    expect_equal(metrics$n_lipids, 3)
    expect_named(
        metrics,
        c(
            "n_lipids",
            "detected_per_sample",
            "missingness",
            "sum_intensity_per_sample",
            "cv_distribution",
            "is_metrics"
        )
    )

    expect_equal(metrics$detected_per_sample$Run, c("S1", "S2", "S3"))
    expect_equal(metrics$detected_per_sample$n_detected, c(2L, 2L, 1L))
    expect_equal(metrics$sum_intensity_per_sample$sum_intensity, c(15, 25, 20))
    expect_equal(metrics$missingness, 4 / 9 * 100)

    l1_group_a <- metrics$cv_distribution[
        metrics$cv_distribution$lipid_id == "L1" &
            metrics$cv_distribution$group == "A",
        ,
        drop = FALSE
    ]
    expect_equal(nrow(l1_group_a), 1)
    expect_equal(l1_group_a$cv, 0)
    expect_equal(nrow(metrics$is_metrics), 0)
})

test_that("calculateLipidFilteringAssayMetrics keeps the non-data-frame placeholder contract", {
    skipIfMissingLipidQcBindings("calculateLipidFilteringAssayMetrics")

    metrics <- NULL
    expect_warning(
        metrics <- calculateLipidFilteringAssayMetrics(
            current_assay_data = 1,
            current_assay_name = "BrokenAssay",
            design_matrix = data.frame(Run = character(), group = character()),
            group_id_col = "group",
            sample_id_col = "Run",
            lipid_id_col = "lipid_id",
            is_pattern = NULL,
            sample_columns = character()
        ),
        "is not a data frame"
    )

    expect_equal(metrics$n_lipids, 0)
    expect_equal(nrow(metrics$detected_per_sample), 0)
    expect_true(is.na(metrics$missingness))
    expect_equal(nrow(metrics$sum_intensity_per_sample), 0)
    expect_equal(nrow(metrics$cv_distribution), 0)
    expect_equal(nrow(metrics$is_metrics), 0)
})

test_that("finalizeLipidFilteringStep updates progress state and returns plot artifacts", {
    skipIfMissingLipidQcBindings(c(
        "calculateLipidFilteringAssayMetrics",
        "finalizeLipidFilteringStep",
        "getFilteringProgressLipidomics"
    ))
    skipIfMissingFilteringProgressLipidomicsClass()

    had_progress <- exists("filtering_progress_lipidomics", envir = .GlobalEnv)
    if (had_progress) {
        original_progress <- get("filtering_progress_lipidomics", envir = .GlobalEnv)
    }
    assign("filtering_progress_lipidomics", new("FilteringProgressLipidomics"), envir = .GlobalEnv)
    on.exit({
        if (had_progress) {
            assign("filtering_progress_lipidomics", original_progress, envir = .GlobalEnv)
        } else if (exists("filtering_progress_lipidomics", envir = .GlobalEnv)) {
            rm("filtering_progress_lipidomics", envir = .GlobalEnv)
        }
    }, add = TRUE)

    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        S1 = c(10, 0, 5),
        S2 = c(10, NA, 15),
        S3 = c(20, 0, 0),
        check.names = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3"),
        group = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )
    metrics <- calculateLipidFilteringAssayMetrics(
        current_assay_data = assay_data,
        current_assay_name = "Assay1",
        design_matrix = design_matrix,
        group_id_col = "group",
        sample_id_col = "Run",
        lipid_id_col = "lipid_id",
        is_pattern = NA_character_,
        sample_columns = c("S1", "S2", "S3")
    )

    step_outputs <- finalizeLipidFilteringStep(
        prog_met = getFilteringProgressLipidomics(),
        step_name = "raw",
        assay_names = "Assay1",
        metrics_list_this_step = list(Assay1 = metrics),
        assay_list = list(Assay1 = assay_data),
        lipid_id_col = "lipid_id",
        overwrite = FALSE
    )

    expect_named(step_outputs, c("total_lipids", "plot_list"))
    expect_equal(step_outputs$total_lipids, 3)
    expect_true(length(step_outputs$plot_list) > 0)
    expect_s3_class(step_outputs$plot_list$total_lipids, "ggplot")

    progress <- get("filtering_progress_lipidomics", envir = .GlobalEnv)
    expect_identical(progress@steps, "raw")
    expect_equal(progress@n_lipids_total, 3)
    expect_named(progress@detected_per_sample[[1]], "Assay1")
})

test_that("saveLipidFilteringPlots persists individual and combined plots", {
    skipIfMissingLipidQcBindings("saveLipidFilteringPlots")

    save_dir <- file.path(tempdir(), paste0("lipid-qc-save-", Sys.getpid()))
    unlink(save_dir, recursive = TRUE, force = TRUE)
    on.exit(unlink(save_dir, recursive = TRUE, force = TRUE), add = TRUE)
    saved_plots <- list()

    plot_list <- list(
        total_lipids = ggplot(mtcars, aes(wt, mpg)) + geom_point(),
        missingness = ggplot(mtcars, aes(factor(cyl), mpg)) + geom_boxplot()
    )

    testthat::local_mocked_bindings(
        ggsave = function(filename, plot, ...) {
            saved_plots[[length(saved_plots) + 1]] <<- list(
                filename = filename,
                plot = plot,
                dots = list(...)
            )
            dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
            file.create(filename)
            invisible(filename)
        },
        .env = environment(saveLipidFilteringPlots)
    )

    saveLipidFilteringPlots(
        plot_list = plot_list,
        step_name = "raw",
        actual_save_dir = save_dir,
        return_grid = TRUE
    )

    expect_true(dir.exists(save_dir))
    expect_true(file.exists(file.path(save_dir, "raw_total_lipids.png")))
    expect_true(file.exists(file.path(save_dir, "raw_missingness.png")))
    expect_true(file.exists(file.path(save_dir, "raw_combined_plots.png")))
    expect_identical(
        vapply(saved_plots, `[[`, character(1), "filename"),
        c(
            file.path(save_dir, "raw_total_lipids.png"),
            file.path(save_dir, "raw_missingness.png"),
            file.path(save_dir, "raw_combined_plots.png")
        )
    )
    expect_s3_class(saved_plots[[3]]$plot, "gtable")
})

test_that("saveLipidFilteringPlots preserves the no-save warning contract", {
    skipIfMissingLipidQcBindings("saveLipidFilteringPlots")

    plot_list <- list(total_lipids = ggplot(mtcars, aes(wt, mpg)) + geom_point())

    expect_warning(
        expect_message(
            saveLipidFilteringPlots(
                plot_list = plot_list,
                step_name = "raw",
                actual_save_dir = NULL,
                return_grid = FALSE,
                publication_graphs_dir = file.path(tempdir(), "requested-output")
            ),
            "Plots will not be saved"
        ),
        "publication_graphs_dir path"
    )
})

test_that("returnLipidFilteringPlots returns a grid grob when requested", {
    skipIfMissingLipidQcBindings("returnLipidFilteringPlots")

    plot_list <- list(
        total_lipids = ggplot(mtcars, aes(wt, mpg)) + geom_point(),
        missingness = ggplot(mtcars, aes(factor(cyl), mpg)) + geom_boxplot()
    )

    grid_plot <- returnLipidFilteringPlots(
        plot_list = plot_list,
        return_grid = TRUE
    )

    expect_false(is.null(grid_plot))
    expect_s3_class(grid_plot, "gtable")
})

test_that("returnLipidFilteringPlots preserves the invisible list-printing contract", {
    skipIfMissingLipidQcBindings("returnLipidFilteringPlots")

    plot_list <- list(
        total_lipids = ggplot(mtcars, aes(wt, mpg)) + geom_point(),
        note = "not-a-plot"
    )

    visible_result <- NULL
    expect_message(
        expect_message(
            visible_result <- withVisible(returnLipidFilteringPlots(
                plot_list = plot_list,
                return_grid = FALSE
            )),
            "Printing plots individually as return_grid is FALSE or grid could not be formed."
        ),
        "Encountered a non-ggplot object"
    )

    expect_false(visible_result$visible)
    expect_identical(visible_result$value, plot_list)
})

test_that("resolveLipidFilteringPlotSaveDir keeps the global project_dirs lookup contract", {
    skipIfMissingLipidQcBindings("resolveLipidFilteringPlotSaveDir")

    tracked_globals <- c("project_dirs", "omic_type", "experiment_label")
    had_globals <- vapply(tracked_globals, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
    original_globals <- lapply(tracked_globals[had_globals], get, envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        for (name in tracked_globals[had_globals]) {
            assign(name, original_globals[[name]], envir = .GlobalEnv)
        }
        for (name in tracked_globals[!had_globals]) {
            if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
                rm(list = name, envir = .GlobalEnv)
            }
        }
    }, add = TRUE)

    assign("project_dirs", list(lipidomics_demo = list(time_dir = tempdir())), envir = .GlobalEnv)
    assign("omic_type", "lipidomics", envir = .GlobalEnv)
    assign("experiment_label", "demo", envir = .GlobalEnv)

    resolved_dir <- NULL
    expect_message(
        resolved_dir <- resolveLipidFilteringPlotSaveDir(
            publication_graphs_dir = file.path(tempdir(), "ignored"),
            omics_type = "ignored_arg",
            time_dir = file.path(tempdir(), "ignored_time_dir")
        ),
        "Attempting to determine save directory"
    )

    expect_identical(resolved_dir, tempdir())
})

test_that("resolveLipidFilteringPlotSaveDir returns NULL when globals are unavailable", {
    skipIfMissingLipidQcBindings("resolveLipidFilteringPlotSaveDir")

    tracked_globals <- c("project_dirs", "omic_type", "experiment_label")
    had_globals <- vapply(tracked_globals, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
    original_globals <- lapply(tracked_globals[had_globals], get, envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        for (name in tracked_globals[had_globals]) {
            assign(name, original_globals[[name]], envir = .GlobalEnv)
        }
        for (name in tracked_globals[!had_globals]) {
            if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
                rm(list = name, envir = .GlobalEnv)
            }
        }
    }, add = TRUE)

    for (name in tracked_globals) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
            rm(list = name, envir = .GlobalEnv)
        }
    }

    expect_warning(
        expect_message(
            expect_null(resolveLipidFilteringPlotSaveDir()),
            "Exists project_dirs: FALSE, Exists omic_type: FALSE, Exists experiment_label: FALSE"
        ),
        "not found"
    )
})

test_that("resolveLipidDuplicateFeaturesByIntensity keeps the highest-average duplicate row", {
    skipIfMissingLipidQcBindings("resolveLipidDuplicateFeaturesByIntensity")

    assay_data <- data.frame(
        lipid_id = c("L1", "L1", "L2"),
        annotation = c("keep", "drop", "solo"),
        S1 = c(10, 1, 5),
        S2 = c(20, 3, 6),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    resolved <- NULL
    expect_message(
        resolved <- resolveLipidDuplicateFeaturesByIntensity(
            assay_tibble = assay_data,
            id_col = "lipid_id",
            sample_cols = c("S1", "S2")
        ),
        "Removed 1 lower-intensity duplicate feature row"
    )

    expect_equal(nrow(resolved), 2)
    expect_equal(resolved$lipid_id, c("L1", "L2"))
    expect_equal(resolved$annotation, c("keep", "solo"))
})

test_that("lipidIntensityFilteringHelper filters rows by sample proportion threshold", {
    skipIfMissingLipidQcBindings("lipidIntensityFilteringHelper")

    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        annotation = c("keep", "drop", "borderline"),
        S1 = c(10, 1, 1),
        S2 = c(12, 2, 9),
        S3 = c(14, 3, 9),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    filtered <- lipidIntensityFilteringHelper(
        assay_table = assay_data,
        min_lipid_intensity_threshold = 5,
        lipids_proportion_of_samples_below_cutoff = 0.5,
        lipid_id_column = "lipid_id"
    )

    expect_equal(filtered$lipid_id, c("L1", "L3"))
    expect_false(any(c("num_below_threshold", "proportion_below_threshold") %in% names(filtered)))
})

test_that("calculateLipidPairCorrelation returns a finite Pearson correlation for one sample pair", {
    skipIfMissingLipidQcBindings("calculateLipidPairCorrelation")

    input_pair_table <- data.frame(
        lipid_id = c("L1", "L1", "L2", "L2", "L3", "L3"),
        Run = c("S1", "S2", "S1", "S2", "S1", "S2"),
        abundance = c(1, 2, 2, 4, 3, 6),
        stringsAsFactors = FALSE
    )

    correlation <- calculateLipidPairCorrelation(
        input_pair_table = input_pair_table,
        feature_id_column = "lipid_id",
        sample_id_column = "Run",
        value_column = "abundance"
    )

    expect_equal(correlation, 1)
})

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

test_that("lipid intensity module output helpers render assay summaries and plots", {
    skipIfMissingLipidQcNamespaceBindings(c(
        "registerLipidIntensityAssayResultsOutput",
        "registerLipidIntensityFilterPlotOutput"
    ))

    output_env <- new.env(parent = emptyenv())
    stats_df <- data.frame(
        Assay = c("AssayA", "AssayB"),
        Original = c(10, 8),
        Filtered = c(8, 5),
        Removed = c(2, 3),
        Percent_Retained = c(80, 62.5),
        stringsAsFactors = FALSE
    )

    lipidQcNamespaceBinding("registerLipidIntensityAssayResultsOutput")(
        output = output_env,
        filterStatsVal = function() stats_df,
        ns = function(id) paste0("intensity-", id),
        renderUiFn = function(expr) htmltools::renderTags(expr)$html
    )

    expect_match(output_env$assay_results_tabs, "intensity-assay_stats_tabs", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "Original lipids:", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "62.5%", fixed = TRUE)

    grob_plot <- structure(list(id = "grid"), class = "grob")
    grob_draws <- list()
    lipidQcNamespaceBinding("registerLipidIntensityFilterPlotOutput")(
        output = output_env,
        filterPlotVal = function() grob_plot,
        renderPlotFn = function(expr) force(expr),
        reqFn = identity,
        gridDrawFn = function(plot_obj) {
            grob_draws <<- c(grob_draws, list(plot_obj))
            invisible(plot_obj)
        },
        printFn = function(...) {
            stop("printFn should not be called for grob plots")
        }
    )
    expect_identical(grob_draws[[1]], grob_plot)

    ggplot_output <- new.env(parent = emptyenv())
    ggplot_plot <- structure(list(id = "gg"), class = "ggplot")
    ggplot_prints <- list()
    lipidQcNamespaceBinding("registerLipidIntensityFilterPlotOutput")(
        output = ggplot_output,
        filterPlotVal = function() ggplot_plot,
        renderPlotFn = function(expr) force(expr),
        reqFn = identity,
        gridDrawFn = function(...) {
            stop("gridDrawFn should not be called for ggplot plots")
        },
        printFn = function(plot_obj) {
            ggplot_prints <<- c(ggplot_prints, list(plot_obj))
            invisible(plot_obj)
        }
    )
    expect_identical(ggplot_prints[[1]], ggplot_plot)
})

test_that("lipid intensity module observers preserve apply and revert contracts", {
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

    lipidQcNamespaceBinding("registerLipidIntensityApplyFilterObserver")(
        input = input,
        workflowData = list(state_manager = state_manager, config_list = list(cache = "cfg")),
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
        lipidIntensityFilteringFn = function(theObject, lipids_intensity_cutoff_percentile, lipids_proportion_of_samples_below_cutoff) {
            filtering_args <<- list(
                theObject = theObject,
                percentile = lipids_intensity_cutoff_percentile,
                proportion = lipids_proportion_of_samples_below_cutoff
            )
            filtered_s4
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

    expect_identical(filtering_args$theObject, current_s4)
    expect_identical(filtering_args$percentile, 10L)
    expect_identical(filtering_args$proportion, 0.5)
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
    expect_identical(save_state_args$state_name, "lipid_intensity_filtered")
    expect_identical(save_state_args$s4_data_object, filtered_s4)
    expect_identical(save_state_args$config_object, list(cache = "cfg"))
    expect_match(output_env$filter_results, "Lipid Intensity Filter Applied Successfully", fixed = TRUE)
    expect_match(output_env$filter_results, "Total: 4 -> 2 lipids (removed 2)", fixed = TRUE)
    expect_true(any(grepl("Lipid intensity filter applied: removed 2 lipids", info_messages, fixed = TRUE)))
    expect_identical(removed_notifications, "lipid_intensity_filter_working")
    expect_identical(notifications[[2]]$message, "Intensity filter applied: 2 lipids retained")

    lipidQcNamespaceBinding("registerLipidIntensityRevertObserver")(
        input = input,
        workflowData = list(state_manager = state_manager),
        output = output_env,
        filterStatsVal = function(value) {
            filter_stats_value <<- value
        },
        filterPlotVal = function(value) {
            filter_plot_value <<- value
        },
        observeEventFn = function(event, handler) {
            expect_identical(event, 8L)
            force(handler)
            invisible(NULL)
        },
        renderTextFn = function(expr) force(expr),
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logErrorFn = function(...) {
            stop("logErrorFn should not be called on revert success")
        },
        showNotificationFn = function(message, type = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type)))
        }
    )

    expect_identical(reverted_state, "lipid_intensity_filtered")
    expect_identical(output_env$filter_results, "Reverted to previous state: lipid_intensity_filtered")
    expect_null(filter_stats_value)
    expect_null(filter_plot_value)
    expect_true(any(grepl("Reverted lipid intensity filter to lipid_intensity_filtered", info_messages, fixed = TRUE)))
})
