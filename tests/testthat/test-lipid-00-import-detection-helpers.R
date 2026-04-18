library(testthat)

test_that("detectLipidomicsFormat distinguishes MS-DIAL and LipidSearch headers", {
    msdial_headers <- c(
        "Peak ID",
        "Lipid name",
        "RT (min)",
        "Precursor m/z",
        "Height",
        "Total score",
        "Sample_A",
        "Sample_B"
    )
    lipidsearch_headers <- c(
        "LipidClass",
        "LipidName",
        "LipidIon",
        "BaseRt",
        "CalcMz",
        "Grade",
        "Sample_A",
        "Sample_B"
    )

    msdial_format <- detectLipidomicsFormat(msdial_headers, filename = "lipid_height_matrix.csv")
    lipidsearch_format <- detectLipidomicsFormat(lipidsearch_headers, filename = "alignment_export.csv")

    expect_identical(msdial_format$format, "msdial")
    expect_true(msdial_format$confidence >= 0.2)
    expect_identical(lipidsearch_format$format, "lipidsearch")
    expect_true(lipidsearch_format$all_scores[["lipidsearch"]] > lipidsearch_format$all_scores[["msdial"]])
})

test_that("getLipidomicsColumnDefaults exposes mirrored vendor defaults", {
    msdial_defaults <- getLipidomicsColumnDefaults("msdial")
    lipidsearch_defaults <- getLipidomicsColumnDefaults("lipidsearch")

    expect_true("Peak ID" %in% msdial_defaults$lipid_id)
    expect_true("Lipid name" %in% msdial_defaults$annotation)
    expect_match(msdial_defaults$is_pattern, "ISTD")

    expect_true("LipidName" %in% lipidsearch_defaults$lipid_id)
    expect_true("LipidClass" %in% lipidsearch_defaults$annotation)
    expect_true("IonType" %in% lipidsearch_defaults$adduct)
})

test_that("findLipidMatchingColumn falls back to the shared matcher", {
    headers <- c("Peak ID", "Lipid name", "Sample_A")

    expect_identical(
        findLipidMatchingColumn(headers, c("peak id", "alignment id")),
        "Peak ID"
    )
    expect_null(findLipidMatchingColumn(headers, c("missing column")))
})

test_that("importLipidMSDIALData keeps numeric annotations out of sample_columns", {
    msdial_file <- tempfile(fileext = ".tsv")
    writeLines(
        c(
            paste(
                c(
                    "Peak ID",
                    "Lipid name",
                    "RT (min)",
                    "Precursor m/z",
                    "Height",
                    "Area",
                    "Total score",
                    "Sample_A",
                    "Sample_B"
                ),
                collapse = "\t"
            ),
            paste(c("1", "PC 34:1", "1.25", "760.6", "100", "120", "88", "101", "99"), collapse = "\t"),
            paste(c("2", "TG 52:2", "2.50", "880.8", "95", "115", "76", "0", "50"), collapse = "\t")
        ),
        con = msdial_file
    )
    on.exit(unlink(msdial_file), add = TRUE)

    import_result <- importLipidMSDIALData(msdial_file)

    expect_identical(import_result$format, "msdial")
    expect_identical(import_result$detected_columns$lipid_id, "Peak ID")
    expect_identical(import_result$detected_columns$annotation, "Lipid name")
    expect_identical(import_result$sample_columns, c("Sample_A", "Sample_B"))
    expect_true(all(c("Height", "Area", "Total score") %in% import_result$annotation_columns))
    expect_match(import_result$is_pattern, "ISTD")
})

test_that("importLipidSearchData keeps metadata columns out of sample_columns", {
    lipidsearch_file <- tempfile(fileext = ".csv")
    writeLines(
        c(
            paste(
                c(
                    "LipidClass",
                    "LipidName",
                    "LipidIon",
                    "BaseRt",
                    "CalcMz",
                    "Grade",
                    "Sample_A",
                    "Sample_B"
                ),
                collapse = ","
            ),
            paste(c("PC", "PC 34:1", "[M+H]+", "1.2", "760.6", "A", "125", "130"), collapse = ","),
            paste(c("TG", "TG 52:2", "[M+NH4]+", "2.8", "880.8", "B", "80", "60"), collapse = ",")
        ),
        con = lipidsearch_file
    )
    on.exit(unlink(lipidsearch_file), add = TRUE)

    import_result <- importLipidSearchData(lipidsearch_file)

    expect_identical(import_result$format, "lipidsearch")
    expect_identical(import_result$detected_columns$lipid_id, "LipidName")
    expect_identical(import_result$detected_columns$annotation, "LipidClass")
    expect_identical(import_result$sample_columns, c("Sample_A", "Sample_B"))
    expect_true(all(c("BaseRt", "CalcMz") %in% import_result$annotation_columns))
})

test_that("validateLipidColumnMapping returns summary stats and missing-column errors", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L2"),
        Sample_A = c(10, 0, 5),
        Sample_B = c(20, NA, 25),
        check.names = FALSE
    )

    valid_mapping <- validateLipidColumnMapping(
        assay_data,
        lipid_id_column = "lipid_id",
        sample_columns = c("Sample_A", "Sample_B")
    )
    invalid_mapping <- validateLipidColumnMapping(
        assay_data,
        lipid_id_column = "missing_id",
        sample_columns = c("Sample_A", "Missing_Sample")
    )

    expect_true(valid_mapping$valid)
    expect_identical(valid_mapping$summary$n_lipids, 2L)
    expect_identical(valid_mapping$summary$n_samples, 2L)
    expect_identical(valid_mapping$summary$pct_missing, 33.3)
    expect_true(any(grepl("duplicate lipid IDs", valid_mapping$warnings)))

    expect_false(invalid_mapping$valid)
    expect_true(any(grepl("Lipid ID column not found", invalid_mapping$errors)))
    expect_true(any(grepl("Sample columns not found", invalid_mapping$errors)))
})
