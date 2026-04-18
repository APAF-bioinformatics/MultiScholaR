library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

test_that("correlation methods resolve through the live lipid S4 correlation helper", {
    wrapper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_objects.R"))
    helper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_correlation_methods.R"))

    expect_false(any(grepl("f = \"pearsonCorForSamplePairs\"", wrapper_source_lines, fixed = TRUE)))
    expect_true(any(grepl("f = \"pearsonCorForSamplePairs\"", helper_source_lines, fixed = TRUE)))
    expect_false(any(grepl("f = \"plotPearson\"", wrapper_source_lines, fixed = TRUE)))
    expect_true(any(grepl("f = \"plotPearson\"", helper_source_lines, fixed = TRUE)))
    expect_false(any(grepl("f = \"filterSamplesByLipidCorrelationThreshold\"", wrapper_source_lines, fixed = TRUE)))
    expect_true(any(grepl("f = \"filterSamplesByLipidCorrelationThreshold\"", helper_source_lines, fixed = TRUE)))
})

test_that("pearsonCorForSamplePairs returns pairwise correlations for lipid technical replicates", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        S1 = c(1, 2, 3),
        S2 = c(2, 4, 6),
        S3 = c(5, 10, 15),
        S4 = c(10, 20, 30),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        tech_rep = c("pair_1", "pair_1", "pair_2", "pair_2"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep"
    )

    correlation_results <- NULL
    expect_warning(
        correlation_results <- pearsonCorForSamplePairs(
            lipid_object,
            tech_rep_remove_regex = ""
        ),
        "many-to-many"
    )

    expect_named(correlation_results, "Assay1")
    expect_s3_class(correlation_results$Assay1, "data.frame")
    expect_identical(correlation_results$Assay1$tech_rep, c("pair_1", "pair_2"))
    expect_identical(correlation_results$Assay1$Run.x, c("S2", "S4"))
    expect_identical(correlation_results$Assay1$Run.y, c("S1", "S3"))
    expect_equal(correlation_results$Assay1$pearson_correlation, c(1, 1))
})

test_that("plotPearson returns named ggplot histograms for lipid technical replicates", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        S1 = c(1, 2, 3),
        S2 = c(2, 4, 6),
        S3 = c(5, 10, 15),
        S4 = c(10, 20, 30),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        tech_rep = c("pair_1", "pair_1", "pair_2", "pair_2"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep"
    )

    pearson_plots <- NULL
    expect_warning(
        pearson_plots <- plotPearson(
            lipid_object,
            tech_rep_remove_regex = ""
        ),
        "many-to-many"
    )

    expect_named(pearson_plots, "Assay1")
    expect_length(pearson_plots, 1)
    expect_s3_class(pearson_plots$Assay1, "ggplot")
    expect_equal(ggplot_build(pearson_plots$Assay1)$plot$labels$x, "Pearson Correlation")
    expect_equal(ggplot_build(pearson_plots$Assay1)$plot$labels$y, "Counts")
})

test_that("filterSamplesByLipidCorrelationThreshold removes low-correlation samples and updates design", {
    assay_data <- data.frame(
        lipid_id = c("L1", "L2", "L3"),
        S1 = c(1, 2, 3),
        S2 = c(2, 4, 6),
        S3 = c(5, 10, 15),
        S4 = c(10, 20, 30),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        tech_rep = c("pair_1", "pair_1", "pair_2", "pair_2"),
        stringsAsFactors = FALSE
    )

    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(Assay1 = assay_data),
        design_matrix = design_matrix,
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep"
    )

    correlation_results <- list(
        Assay1 = data.frame(
            Run.x = c("S2", "S4"),
            Run.y = c("S1", "S3"),
            pearson_correlation = c(0.9, 0.2),
            stringsAsFactors = FALSE
        )
    )

    filtered_object <- filterSamplesByLipidCorrelationThreshold(
        lipid_object,
        pearson_correlation_per_pair = correlation_results,
        min_pearson_correlation_threshold = 0.5
    )

    expect_identical(colnames(filtered_object@lipid_data$Assay1), c("lipid_id", "S1", "S2"))
    expect_identical(filtered_object@design_matrix$Run, c("S1", "S2"))
})
