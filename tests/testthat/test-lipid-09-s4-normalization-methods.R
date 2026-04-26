# fidelity-coverage-compare: shared
library(testthat)
suppressPackageStartupMessages(library(dplyr))

assign("log_info", logger::log_info, envir = globalenv())
assign("log_warn", logger::log_warn, envir = globalenv())

expect_method_defined_in_current_layout <- function(method_name, wrapper_path, helper_path) {
    wrapper_source_lines <- readLines(wrapper_path)
    wrapper_has_method <- any(grepl(sprintf('f = "%s"', method_name), wrapper_source_lines, fixed = TRUE))
    helper_has_method <- FALSE

    if (file.exists(helper_path)) {
        helper_source_lines <- readLines(helper_path)
        helper_has_method <- any(grepl(sprintf('f = "%s"', method_name), helper_source_lines, fixed = TRUE))
    }

    expect_true(wrapper_has_method || helper_has_method)

    if (file.exists(helper_path)) {
        expect_false(wrapper_has_method)
        expect_true(helper_has_method)
    } else {
        expect_true(wrapper_has_method)
    }
}

test_that("normaliseUntransformedData records ITSD normalization outputs", {
    assay_data <- tibble::tibble(
        Name = c("ITSD_1", "L1", "L2"),
        Annotation = c("ITSD_1", "Alpha", "Beta"),
        S1 = c(10, 100, 40),
        S2 = c(20, 100, 80)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group",
        internal_standard_regex = "^ITSD"
    )

    normalized_object <- normaliseUntransformedData(
        lipid_object,
        method = "ITSD",
        remove_itsd_after_norm = TRUE
    )

    normalized_assay <- normalized_object@lipid_data$Assay1

    expect_equal(normalized_assay$Name, c("L1", "L2"))
    expect_equal(normalized_assay$S1, c(150, 60))
    expect_equal(normalized_assay$S2, c(75, 60))
    expect_true(isTRUE(normalized_object@args$ITSDNormalization$applied))
    expect_identical(normalized_object@args$ITSDNormalization$itsd_aggregation, "sum")
    expect_identical(normalized_object@args$ITSDNormalization$itsd_pattern_columns, "Annotation")
    expect_equal(normalized_object@args$ITSDNormalization$itsd_counts_per_assay$Assay1, 1)
    expect_identical(normalized_object@args$ITSDNormalization$itsd_features_per_assay$Assay1, "ITSD_1")
})

test_that("lipidIntensityFiltering resolves through the live lipid S4 normalization helper", {
    expect_method_defined_in_current_layout(
        method_name = "lipidIntensityFiltering",
        wrapper_path = file.path("..", "..", "R", "func_lipid_s4_objects.R"),
        helper_path = file.path("..", "..", "R", "func_lipid_s4_normalization_methods.R")
    )
})

test_that("lipidIntensityFiltering filters low-intensity lipids across assays", {
    assay_data <- tibble::tibble(
        Name = c("L1", "L2", "L3"),
        Annotation = c("Alpha", "Beta", "Gamma"),
        S1 = c(1, 10, 100),
        S2 = c(1, 10, 100)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    filtered_object <- lipidIntensityFiltering(
        lipid_object,
        lipids_intensity_cutoff_percentile = 50,
        lipids_proportion_of_samples_below_cutoff = 0.5
    )

    expect_identical(filtered_object@lipid_data$Assay1$Name, c("L2", "L3"))
    expect_identical(colnames(filtered_object@lipid_data$Assay1), colnames(assay_data))
})

test_that("getNegCtrlMetabAnova resolves through the live lipid RUV helper and short-circuits empty assays", {
    expect_method_defined_in_current_layout(
        method_name = "getNegCtrlMetabAnova",
        wrapper_path = file.path("..", "..", "R", "func_lipid_s4_objects.R"),
        helper_path = file.path("..", "..", "R", "func_lipid_norm_ruv_helpers.R")
    )

    empty_object <- createLipidomicsAssayData(
        lipid_data = list(),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    expect_identical(getNegCtrlMetabAnova(empty_object), list())
})

test_that("ruvCancor resolves through the live lipid RUV helper and rejects missing controls early", {
    expect_method_defined_in_current_layout(
        method_name = "ruvCancor",
        wrapper_path = file.path("..", "..", "R", "func_lipid_s4_objects.R"),
        helper_path = file.path("..", "..", "R", "func_lipid_norm_ruv_helpers.R")
    )

    empty_object <- createLipidomicsAssayData(
        lipid_data = list(),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    expect_error(
        ruvCancor(empty_object),
        "ctrl' is not defined|Missing required 'ctrl' parameter for ruvCancor"
    )
})

test_that("ruvIII_C_Varying resolves through the live lipid RUV helper and rejects missing k early", {
    expect_method_defined_in_current_layout(
        method_name = "ruvIII_C_Varying",
        wrapper_path = file.path("..", "..", "R", "func_lipid_s4_objects.R"),
        helper_path = file.path("..", "..", "R", "func_lipid_norm_ruv_helpers.R")
    )

    empty_object <- createLipidomicsAssayData(
        lipid_data = list(),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    expect_error(
        ruvIII_C_Varying(empty_object, ruv_grouping_variable = "group"),
        "ruv_number_k' is not defined|Missing required 'ruv_number_k'"
    )
})

test_that("normaliseBetweenSamples preserves assay values for none and cleans design order", {
    assay_data <- tibble::tibble(
        Name = c("L1", "L2"),
        Annotation = c("Alpha", "Beta"),
        S2 = c(20, 40),
        S1 = c(10, 30)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )
    lipid_object@args$normaliseBetweenSamples <- list(normalisation_method = "none")

    normalized_object <- normaliseBetweenSamples(lipid_object)

    expect_equal(normalized_object@lipid_data$Assay1, assay_data)
    expect_identical(normalized_object@args$normalisation_method, "none")
    expect_identical(normalized_object@design_matrix$Run, c("S2", "S1"))
    expect_identical(normalized_object@design_matrix$group, c("B", "A"))
})

test_that("cleanDesignMatrix drops unmatched design rows and follows assay order", {
    assay_data <- tibble::tibble(
        Name = c("L1", "L2"),
        Annotation = c("Alpha", "Beta"),
        S2 = c(20, 40),
        S1 = c(10, 30)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        batch = c("b1", "b2"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )
    lipid_object@design_matrix <- data.frame(
        Run = c("S3", "S1", "S2"),
        group = c("C", "A", "B"),
        batch = c("b3", "b1", "b2"),
        stringsAsFactors = FALSE
    )

    cleaned_object <- cleanDesignMatrix(lipid_object)

    expect_identical(cleaned_object@design_matrix$Run, c("S2", "S1"))
    expect_identical(cleaned_object@design_matrix$group, c("B", "A"))
    expect_identical(cleaned_object@design_matrix$batch, c("b2", "b1"))
})

test_that("logTransformAssays log2 transforms sample columns and records args", {
    assay_data <- tibble::tibble(
        Name = c("L1", "L2"),
        Annotation = c("Alpha", "Beta"),
        Class = c("C1", "C2"),
        S1 = c(3, -2),
        S2 = c(0, 7)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2"),
        group = c("A", "B"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    transformed_object <- logTransformAssays(lipid_object, offset = 1)
    transformed_assay <- transformed_object@lipid_data$Assay1

    expect_identical(transformed_assay$Name, assay_data$Name)
    expect_identical(transformed_assay$Annotation, assay_data$Annotation)
    expect_identical(transformed_assay$Class, assay_data$Class)
    expect_equal(transformed_assay$S1, c(log2(4), log2(1)))
    expect_equal(transformed_assay$S2, c(log2(1), log2(8)))
    expect_true(isTRUE(transformed_object@args$log_transformed))
    expect_identical(transformed_object@args$log_transform_offset, 1)
})
