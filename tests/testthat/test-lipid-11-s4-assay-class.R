# fidelity-coverage-compare: shared
library(testthat)

test_that("LipidomicsAssayData exposes the expected slot contract", {
    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L2"),
                S1 = c(10, 20),
                S2 = c(11, 21),
                check.names = FALSE
            )
        ),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group"
    )

    expect_s4_class(lipid_object, "LipidomicsAssayData")
    expect_identical(methods::slotNames(lipid_object), c(
        "lipid_data",
        "lipid_id_column",
        "annotation_id_column",
        "database_identifier_type",
        "internal_standard_regex",
        "design_matrix",
        "sample_id",
        "group_id",
        "technical_replicate_id",
        "args"
    ))
    expect_true(isTRUE(methods::validObject(lipid_object, test = TRUE)))
})

test_that("LipidomicsAssayData validity keeps the missing-sample-id-column error", {
    expect_error(
        methods::new(
            "LipidomicsAssayData",
            lipid_data = list(
                AssayA = data.frame(
                    lipid_id = c("L1", "L2"),
                    S1 = c(10, 20),
                    S2 = c(11, 21),
                    check.names = FALSE
                )
            ),
            lipid_id_column = "lipid_id",
            design_matrix = data.frame(
                Sample = c("S1", "S2"),
                group = c("A", "B"),
                stringsAsFactors = FALSE
            ),
            sample_id = "Run",
            group_id = "group"
        ),
        regexp = "`sample_id` column \\('Run'\\) not found"
    )
})

test_that("LipidomicsAssayData validity keeps assay/design mismatch details", {
    expect_error(
        methods::new(
            "LipidomicsAssayData",
            lipid_data = list(
                AssayA = data.frame(
                    lipid_id = c("L1", "L2"),
                    S1 = c(10, 20),
                    S3 = c(11, 21),
                    check.names = FALSE
                )
            ),
            lipid_id_column = "lipid_id",
            design_matrix = data.frame(
                Run = c("S1", "S2"),
                group = c("A", "B"),
                stringsAsFactors = FALSE
            ),
            sample_id = "Run",
            group_id = "group"
        ),
        regexp = "Sample columns in assays do not exactly match unique sample IDs \\('Run'\\)"
    )
    expect_error(
        methods::new(
            "LipidomicsAssayData",
            lipid_data = list(
                AssayA = data.frame(
                    lipid_id = c("L1", "L2"),
                    S1 = c(10, 20),
                    S3 = c(11, 21),
                    check.names = FALSE
                )
            ),
            lipid_id_column = "lipid_id",
            design_matrix = data.frame(
                Run = c("S1", "S2"),
                group = c("A", "B"),
                stringsAsFactors = FALSE
            ),
            sample_id = "Run",
            group_id = "group"
        ),
        regexp = "Samples in design_matrix missing from assay columns: S2"
    )
})
