library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

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

    had_runner <- exists("runTestsContrasts", envir = .GlobalEnv, inherits = FALSE)
    old_runner <- if (had_runner) {
        get("runTestsContrasts", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had_runner) {
            assign("runTestsContrasts", old_runner, envir = .GlobalEnv)
        } else if (exists("runTestsContrasts", envir = .GlobalEnv, inherits = FALSE)) {
            rm("runTestsContrasts", envir = .GlobalEnv)
        }
    }, add = TRUE)

    captured_call <- NULL
    assign(
        "runTestsContrasts",
        function(data_matrix, contrast_strings, design_matrix, formula_string, treat_lfc_cutoff, eBayes_trend, eBayes_robust) {
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
        },
        envir = .GlobalEnv
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
