library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(stringr))

test_that("getDaResults accessor methods populate lipid DA result tables", {
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

    wide_results <- getDaResultsWideFormat(list(Assay1 = da_results_object))
    expect_named(wide_results, "Assay1")
    expect_s4_class(wide_results$Assay1, "LipidomicsDifferentialAbundanceResults")
    expect_true(all(c("lipid_id", "annotation", "S1", "S2") %in% colnames(wide_results$Assay1@results_table_wide)))
    expect_equal(nrow(wide_results$Assay1@results_table_wide), 2)
    expect_equal(
        wide_results$Assay1@results_table_wide$S1[
            wide_results$Assay1@results_table_wide$lipid_id == "L1"
        ],
        10
    )
    expect_true(any(grepl("^logFC:", colnames(wide_results$Assay1@results_table_wide))))

    long_results <- getDaResultsLongFormat(list(Assay1 = da_results_object))
    expect_named(long_results, "Assay1")
    expect_s4_class(long_results$Assay1, "LipidomicsDifferentialAbundanceResults")
    expect_true(all(c("comparison", "lipid_id", "logFC", "annotation", "S1", "S2") %in% colnames(long_results$Assay1@results_table_long)))
    expect_equal(nrow(long_results$Assay1@results_table_long), 2)
    expect_equal(
        long_results$Assay1@results_table_long$logFC[
            long_results$Assay1@results_table_long$lipid_id == "L2"
        ],
        -1.25
    )
})
