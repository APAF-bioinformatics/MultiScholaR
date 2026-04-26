# fidelity-coverage-compare: shared
library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

localLipidDaPlotBinding <- function(env, name, value, .local_envir = parent.frame()) {
    had_binding <- exists(name, envir = env, inherits = FALSE)
    old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
    was_locked <- had_binding && bindingIsLocked(name, env)

    if (was_locked) {
        unlockBinding(name, env)
    }
    assign(name, value, envir = env)
    if (was_locked) {
        lockBinding(name, env)
    }

    withr::defer({
        if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
            unlockBinding(name, env)
        }
        if (had_binding) {
            assign(name, old_value, envir = env)
        } else if (exists(name, envir = env, inherits = FALSE)) {
            rm(list = name, envir = env)
        }
        if (was_locked && exists(name, envir = env, inherits = FALSE)) {
            lockBinding(name, env)
        }
    }, envir = .local_envir)
}

test_that("plotNumSigDiffExpBarPlot stores the summary plot and table on each assay", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2"),
        annotation = c("Alpha", "Beta"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group"
    )

    da_results_object <- new(
        "LipidomicsDifferentialAbundanceResults",
        theObject = lipid_object,
        contrasts_results_table = list(
            "groupB-groupA" = data.frame(
                lipid_id = c("L1", "L2"),
                logFC = c(1.5, -1.25),
                raw_pvalue = c(0.01, 0.02),
                fdr_qvalue = c(0.015, 0.025),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        )
    )

    captured_descriptions <- NULL
    counter_mock <- function(list_of_da_tables, list_of_descriptions, formula_string) {
        captured_descriptions <<- list_of_descriptions
        expect_type(list_of_da_tables, "list")
        expect_identical(formula_string, NA)
        list(
            plot = list(summary_plot = "mock-bar-plot"),
            table = data.frame(
                comparison = list_of_descriptions,
                significant = 2L,
                stringsAsFactors = FALSE
            )
        )
    }
    localLipidDaPlotBinding(
        .GlobalEnv,
        "printCountDaGenesTable",
        counter_mock
    )
    localLipidDaPlotBinding(
        asNamespace("MultiScholaR"),
        "printCountDaGenesTable",
        counter_mock
    )

    plotted_results <- plotNumSigDiffExpBarPlot(list(Assay1 = da_results_object))

    expect_named(plotted_results, "Assay1")
    expect_identical(captured_descriptions, "groupB-groupA")
    expect_identical(plotted_results$Assay1@num_sig_diff_exp_bar_plot$summary_plot, "mock-bar-plot")
    expect_identical(plotted_results$Assay1@num_sig_diff_table$comparison, "groupB-groupA")
    expect_identical(plotted_results$Assay1@num_sig_diff_table$significant, 2L)
})

test_that("plotVolcanoS4 stores one named plot per contrast with classified point labels", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2"),
        annotation = c("Alpha", "Beta"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group"
    )

    da_results_object <- new(
        "LipidomicsDifferentialAbundanceResults",
        theObject = lipid_object,
        contrasts_results_table = list(
            "groupB-groupA" = data.frame(
                lipid_id = c("L1", "L2"),
                logFC = c(1.5, -1.25),
                raw_pvalue = c(0.01, 0.02),
                fdr_qvalue = c(0.015, 0.025),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        )
    )

    volcano_calls <- list()
    volcano_mock <- function(input_data, input_title, log_q_value_column, log_fc_column) {
        volcano_calls <<- append(volcano_calls, list(list(
            data = input_data,
            title = input_title,
            log_q_value_column = log_q_value_column,
            log_fc_column = log_fc_column
        )))
        ggplot2::ggplot(input_data, ggplot2::aes(.data[[log_fc_column]], .data[[log_q_value_column]], colour = label)) +
            ggplot2::geom_point()
    }
    localLipidDaPlotBinding(
        .GlobalEnv,
        "plotOneVolcanoNoVerticalLines",
        volcano_mock
    )
    localLipidDaPlotBinding(
        asNamespace("MultiScholaR"),
        "plotOneVolcanoNoVerticalLines",
        volcano_mock
    )

    plotted_results <- plotVolcanoS4(list(Assay1 = da_results_object))

    expect_named(plotted_results, "Assay1")
    expect_named(plotted_results$Assay1@volcano_plot, "groupB-groupA")
    expect_length(volcano_calls, 1)
    expect_identical(volcano_calls[[1]]$title, "Assay1 - groupB-groupA")
    expect_identical(volcano_calls[[1]]$log_q_value_column, "lqm")
    expect_identical(volcano_calls[[1]]$log_fc_column, "logFC")
    expect_identical(as.character(volcano_calls[[1]]$data$label), c("Significant Up", "Significant Down"))
    expect_identical(as.character(volcano_calls[[1]]$data$colour), c("red", "blue"))
    expect_s3_class(plotted_results$Assay1@volcano_plot[[1]], "ggplot")
})
