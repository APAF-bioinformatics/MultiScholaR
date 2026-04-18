library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

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

    had_counter <- exists("printCountDaGenesTable", envir = .GlobalEnv, inherits = FALSE)
    old_counter <- if (had_counter) {
        get("printCountDaGenesTable", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had_counter) {
            assign("printCountDaGenesTable", old_counter, envir = .GlobalEnv)
        } else if (exists("printCountDaGenesTable", envir = .GlobalEnv, inherits = FALSE)) {
            rm("printCountDaGenesTable", envir = .GlobalEnv)
        }
    }, add = TRUE)

    captured_descriptions <- NULL
    assign(
        "printCountDaGenesTable",
        function(list_of_da_tables, list_of_descriptions, formula_string) {
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
        },
        envir = .GlobalEnv
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

    restore_name <- function(name, old_value, existed) {
        if (existed) {
            assign(name, old_value, envir = .GlobalEnv)
        } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
            rm(list = name, envir = .GlobalEnv)
        }
    }

    had_volcano <- exists("plotOneVolcanoNoVerticalLines", envir = .GlobalEnv, inherits = FALSE)
    old_volcano <- if (had_volcano) {
        get("plotOneVolcanoNoVerticalLines", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    had_scale <- exists("scale_color_manual", envir = .GlobalEnv, inherits = FALSE)
    old_scale <- if (had_scale) {
        get("scale_color_manual", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit(restore_name("plotOneVolcanoNoVerticalLines", old_volcano, had_volcano), add = TRUE)
    on.exit(restore_name("scale_color_manual", old_scale, had_scale), add = TRUE)

    volcano_calls <- list()
    assign(
        "plotOneVolcanoNoVerticalLines",
        function(input_data, input_title, log_q_value_column, log_fc_column) {
            volcano_calls <<- append(volcano_calls, list(list(
                data = input_data,
                title = input_title,
                log_q_value_column = log_q_value_column,
                log_fc_column = log_fc_column
            )))
            ggplot2::ggplot()
        },
        envir = .GlobalEnv
    )
    assign(
        "scale_color_manual",
        function(...) {
            NULL
        },
        envir = .GlobalEnv
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
