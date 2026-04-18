library(testthat)

test_that("getCountsTable resolves lipid assay counts through the canonical DA helper", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        check.names = FALSE
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
        sample_id = "Run",
        group_id = "group"
    )

    expect_identical(getCountsTable(lipid_object), lipid_object@lipid_data)
    expect_error(getCountsTable(list()), "Unsupported object type")
})

test_that("runTestsContrastsLipidDA rejects undefined contrast levels early", {
    expr_matrix <- matrix(
        c(
            10, 12, 30, 32,
            20, 18, 5, 7
        ),
        nrow = 2,
        byrow = TRUE
    )
    rownames(expr_matrix) <- c("L1", "L2")
    colnames(expr_matrix) <- c("S1", "S2", "S3", "S4")

    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    expect_error(
        runTestsContrastsLipidDA(
            data = expr_matrix,
            contrast_strings = "Treatment-Control",
            design_matrix = design_matrix,
            formula_string = "~ 0 + group",
            sample_id_col = "Run",
            treat_lfc_cutoff = NA
        ),
        "references undefined levels"
    )
})

test_that("createLipidDaResultsLongFormat appends intensities and parsed groups", {
    expr_matrix <- matrix(
        c(
            10, 12, 30, 32,
            20, 18, 5, 7
        ),
        nrow = 2,
        byrow = TRUE
    )
    rownames(expr_matrix) <- c("L1", "L2")
    colnames(expr_matrix) <- c("S1", "S2", "S3", "S4")

    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    long_tbl <- createLipidDaResultsLongFormat(
        lfc_qval_tbl = data.frame(
            lipid_id = c("L1", "L2"),
            comparison = c("groupB-groupA", "groupB-groupA"),
            logFC = c(1.0, -1.0),
            raw_pvalue = c(0.01, 0.02),
            fdr_qvalue = c(0.02, 0.03),
            stringsAsFactors = FALSE
        ),
        expr_matrix = expr_matrix,
        design_matrix = design_matrix,
        sample_id_col = "Run",
        group_id_col = "group",
        lipid_id_col = "lipid_id"
    )

    expect_equal(nrow(long_tbl), 2)
    expect_true(all(c("intensity.S1.A", "intensity.S4.B") %in% colnames(long_tbl)))
    expect_identical(unique(long_tbl$numerator), "groupB")
    expect_identical(unique(long_tbl$denominator), "groupA")
    expect_equal(long_tbl$intensity.S1.A[long_tbl$lipid_id == "L1"], 10)
    expect_equal(long_tbl$intensity.S4.B[long_tbl$lipid_id == "L2"], 7)
})

test_that("runLipidsDA returns long-format assay results with intensity columns", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2"),
        lipid_name = c("Lipid 1", "Lipid 2"),
        S1 = c(10, 20),
        S2 = c(12, 18),
        S3 = c(30, 5),
        S4 = c(32, 7),
        check.names = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        annotation_id_column = "lipid_name",
        sample_id = "Run",
        group_id = "group"
    )

    original_runner <- runTestsContrastsLipidDA
    assign(
        "runTestsContrastsLipidDA",
        function(...) {
            list(
                results = list(
                    "groupB-groupA" = data.frame(
                        logFC = c(1.8, -1.6),
                        P.Value = c(0.01, 0.02),
                        fdr_qvalue = c(0.02, 0.03),
                        raw_pvalue = c(0.01, 0.02),
                        row.names = c("L1", "L2"),
                        check.names = FALSE
                    )
                ),
                fit.eb = list(coefficients = matrix(c(1.8, -1.6), ncol = 1)),
                qvalue_warnings = list("groupB-groupA" = TRUE)
            )
        },
        envir = .GlobalEnv
    )
    on.exit(
        assign("runTestsContrastsLipidDA", original_runner, envir = .GlobalEnv),
        add = TRUE
    )

    results <- runLipidsDA(
        theObject = lipid_object,
        contrasts_tbl = data.frame(
            contrasts = "groupB-groupA",
            friendly_names = "B_vs_A",
            stringsAsFactors = FALSE
        ),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0
    )

    expect_named(
        results,
        c(
            "theObject",
            "contrasts_results",
            "da_lipids_long",
            "per_assay_results",
            "significant_counts",
            "qvalue_warnings",
            "da_q_val_thresh",
            "treat_lfc_cutoff"
        )
    )
    expect_true("Assay1" %in% names(results$per_assay_results))
    expect_equal(nrow(results$da_lipids_long), 2)
    expect_true(all(c("intensity.S1.A", "intensity.S4.B") %in% colnames(results$da_lipids_long)))
    expect_identical(unique(results$da_lipids_long$friendly_name), "B_vs_A")
    expect_identical(unique(results$da_lipids_long$assay), "Assay1")
    expect_identical(names(results$qvalue_warnings$Assay1), "groupB-groupA")
})
