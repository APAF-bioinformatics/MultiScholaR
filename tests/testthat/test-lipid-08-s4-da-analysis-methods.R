# fidelity-coverage-compare: shared
library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

localLipidDaAnalysisBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

test_that("lipid S4 differentialAbundanceAnalysis methods return named DA result objects", {
    assay_data <- data.frame(
        Name = c("L1", "L2"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        S3 = c(30, 5),
        S4 = c(32, 7),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        sample_id = "Run",
        group_id = "group"
    )

    captured_call <- NULL
    run_tests_mock <- function(data_matrix, contrast_strings, design_matrix, formula_string, treat_lfc_cutoff, eBayes_trend, eBayes_robust) {
        captured_call <<- list(
            rownames = rownames(data_matrix),
            colnames = colnames(data_matrix),
            contrast_strings = contrast_strings,
            formula_string = formula_string,
            treat_lfc_cutoff = treat_lfc_cutoff,
            eBayes_trend = eBayes_trend,
            eBayes_robust = eBayes_robust,
            design_groups = design_matrix$group
        )

        list(
            fit.eb = list(coefficients = matrix(c(1.8, -1.6), ncol = 1)),
            results = list(
                "groupB-groupA" = data.frame(
                    logFC = c(1.8, -1.6),
                    raw_pvalue = c(0.01, 0.02),
                    fdr_qvalue = c(0.02, 0.03),
                    row.names = c("L1", "L2"),
                    check.names = FALSE
                )
            )
        )
    }
    localLipidDaAnalysisBinding(
        .GlobalEnv,
        "runTestsContrasts",
        run_tests_mock
    )
    localLipidDaAnalysisBinding(
        asNamespace("MultiScholaR"),
        "runTestsContrasts",
        run_tests_mock
    )

    helper_result <- differentialAbundanceAnalysisHelper(
        lipid_object,
        contrasts_tbl = data.frame(
            contrasts = "groupB-groupA",
            stringsAsFactors = FALSE
        ),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0
    )

    expect_s4_class(helper_result, "LipidomicsDifferentialAbundanceResults")
    expect_equal(helper_result@contrasts_results_table[[1]]$Name, c("L1", "L2"))
    expect_identical(captured_call$rownames, c("L1", "L2"))
    expect_identical(captured_call$colnames, c("S1", "S2", "S3", "S4"))
    expect_identical(captured_call$contrast_strings, "groupB-groupA")
    expect_identical(captured_call$design_groups, c("A", "A", "B", "B"))

    wrapper_result <- differentialAbundanceAnalysis(
        list(Assay1 = lipid_object),
        contrasts_tbl = data.frame(
            contrasts = "groupB-groupA",
            stringsAsFactors = FALSE
        ),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0
    )

    expect_named(wrapper_result, "Assay1")
    expect_s4_class(wrapper_result$Assay1, "LipidomicsDifferentialAbundanceResults")
    expect_equal(wrapper_result$Assay1@contrasts_results_table[[1]]$Name, c("L1", "L2"))
})

test_that("lipid S4 differentialAbundanceAnalysisHelper sanitizes numeric group names and restores result labels", {
    assay_data <- data.frame(
        Name = c("L1", "L2"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        S3 = c(30, 5),
        S4 = c(32, 7),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("1A", "1A", "2B", "2B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        sample_id = "Run",
        group_id = "group"
    )

    captured_call <- NULL
    run_tests_mock <- function(data_matrix, contrast_strings, design_matrix, formula_string, treat_lfc_cutoff, eBayes_trend, eBayes_robust) {
        captured_call <<- list(
            contrast_strings = contrast_strings,
            formula_string = formula_string,
            treat_lfc_cutoff = treat_lfc_cutoff,
            eBayes_trend = eBayes_trend,
            eBayes_robust = eBayes_robust,
            design_groups = design_matrix$group
        )

        list(
            fit.eb = list(coefficients = matrix(c(1.8, -1.6), ncol = 1)),
            results = list(
                "grp_2B-grp_1A" = data.frame(
                    comparison = c("grp_2B-grp_1A", "grp_2B-grp_1A"),
                    logFC = c(1.8, -1.6),
                    raw_pvalue = c(0.01, 0.02),
                    fdr_qvalue = c(0.02, 0.03),
                    row.names = c("L1", "L2"),
                    check.names = FALSE
                )
            )
        )
    }
    localLipidDaAnalysisBinding(
        .GlobalEnv,
        "runTestsContrasts",
        run_tests_mock
    )
    localLipidDaAnalysisBinding(
        asNamespace("MultiScholaR"),
        "runTestsContrasts",
        run_tests_mock
    )

    helper_result <- differentialAbundanceAnalysisHelper(
        lipid_object,
        contrasts_tbl = data.frame(
            contrasts = "2B-1A",
            stringsAsFactors = FALSE
        ),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0.5,
        eBayes_trend = "FALSE",
        eBayes_robust = "TRUE"
    )

    expect_s4_class(helper_result, "LipidomicsDifferentialAbundanceResults")
    expect_identical(captured_call$contrast_strings, "grp_2B-grp_1A")
    expect_identical(captured_call$design_groups, c("grp_1A", "grp_1A", "grp_2B", "grp_2B"))
    expect_identical(captured_call$treat_lfc_cutoff, 0.5)
    expect_false(captured_call$eBayes_trend)
    expect_true(captured_call$eBayes_robust)
    expect_identical(helper_result@contrasts_results_table[[1]]$comparison, c("2B-1A", "2B-1A"))
    expect_identical(helper_result@theObject@design_matrix$group, c("grp_1A", "grp_1A", "grp_2B", "grp_2B"))
})

test_that("lipid S4 differentialAbundanceAnalysisHelper falls back to namespace runTestsContrasts", {
    assay_data <- data.frame(
        Name = c("L1", "L2"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        S3 = c(30, 5),
        S4 = c(32, 7),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        sample_id = "Run",
        group_id = "group"
    )

    captured_call <- NULL
    run_tests_mock <- function(data_matrix, contrast_strings, design_matrix, formula_string, treat_lfc_cutoff, eBayes_trend, eBayes_robust) {
        captured_call <<- list(
            rownames = rownames(data_matrix),
            contrast_strings = contrast_strings
        )

        list(
            fit.eb = list(coefficients = matrix(c(1.8, -1.6), ncol = 1)),
            results = list(
                "groupB-groupA" = data.frame(
                    logFC = c(1.8, -1.6),
                    raw_pvalue = c(0.01, 0.02),
                    fdr_qvalue = c(0.02, 0.03),
                    row.names = c("L1", "L2"),
                    check.names = FALSE
                )
            )
        )
    }
    localLipidDaAnalysisBinding(
        .GlobalEnv,
        "runTestsContrasts",
        "not-a-function"
    )
    localLipidDaAnalysisBinding(
        asNamespace("MultiScholaR"),
        "runTestsContrasts",
        run_tests_mock
    )

    helper_result <- differentialAbundanceAnalysisHelper(
        lipid_object,
        contrasts_tbl = data.frame(
            contrasts = "groupB-groupA",
            stringsAsFactors = FALSE
        ),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0
    )

    expect_s4_class(helper_result, "LipidomicsDifferentialAbundanceResults")
    expect_identical(captured_call$rownames, c("L1", "L2"))
    expect_identical(captured_call$contrast_strings, "groupB-groupA")
})
