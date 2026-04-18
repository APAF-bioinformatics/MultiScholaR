library(testthat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tibble))

test_that("plotPca resolves through the active lipid S4 plotting source", {
    wrapper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_objects.R"))
    helper_path <- file.path("..", "..", "R", "func_lipid_s4_plotting_methods.R")
    helper_exists <- file.exists(helper_path)

    if (helper_exists) {
        helper_source_lines <- readLines(helper_path)

        expect_false(any(grepl("f = \"plotPca\"", wrapper_source_lines, fixed = TRUE)))
        expect_true(any(grepl("f = \"plotPca\"", helper_source_lines, fixed = TRUE)))
    } else {
        expect_true(any(grepl("f = \"plotPca\"", wrapper_source_lines, fixed = TRUE)))
    }
})

test_that("plotPca returns one titled ggplot per assay with PCA coordinates", {
    assay_data <- tibble::tibble(
        Name = paste0("L", seq_len(10)),
        Annotation = paste0("A", seq_len(10)),
        S1 = c(10, 13, 18, 25, 40, 55, 70, 85, 95, 110),
        S2 = c(11, 14, 17, 24, 42, 57, 69, 84, 93, 108),
        S3 = c(30, 28, 19, 14, 60, 48, 39, 28, 19, 12),
        S4 = c(31, 27, 20, 13, 58, 46, 41, 30, 21, 14)
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
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    had_plot_helper <- exists("plotPcaHelper", envir = .GlobalEnv, inherits = FALSE)
    old_plot_helper <- if (had_plot_helper) {
        get("plotPcaHelper", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had_plot_helper) {
            assign("plotPcaHelper", old_plot_helper, envir = .GlobalEnv)
        } else if (exists("plotPcaHelper", envir = .GlobalEnv, inherits = FALSE)) {
            rm("plotPcaHelper", envir = .GlobalEnv)
        }
    }, add = TRUE)

    captured_call <- NULL
    assign(
        "plotPcaHelper",
        function(data,
                 design_matrix,
                 sample_id_column,
                 grouping_variable,
                 shape_variable = NULL,
                 label_column = NULL,
                 title,
                 geom.text.size = 11,
                 ...) {
            captured_call <<- list(
                data = data,
                design_matrix = design_matrix,
                sample_id_column = sample_id_column,
                grouping_variable = grouping_variable,
                shape_variable = shape_variable,
                label_column = label_column,
                title = title,
                geom.text.size = geom.text.size
            )

            ggplot2::ggplot(
                data.frame(
                    PC1 = c(1, -1),
                    PC2 = c(2, -2),
                    Run = c("S1", "S2"),
                    group = c("A", "B")
                ),
                ggplot2::aes(PC1, PC2, colour = group)
            ) +
                ggplot2::labs(title = title)
        },
        envir = .GlobalEnv
    )

    pca_plots <- plotPca(
        lipid_object,
        grouping_variable = "group",
        title = "QC PCA",
        font_size = 6
    )

    expect_named(pca_plots, "Assay1")
    expect_s3_class(pca_plots$Assay1, "ggplot")
    expect_identical(pca_plots$Assay1$labels$title, "QC PCA - Assay1")
    expect_true(all(c("PC1", "PC2", "Run", "group") %in% colnames(pca_plots$Assay1$data)))
    expect_equal(sort(unique(as.character(pca_plots$Assay1$data$group))), c("A", "B"))
    expect_identical(captured_call$sample_id_column, "Run")
    expect_identical(captured_call$grouping_variable, "group")
    expect_identical(captured_call$title, "QC PCA - Assay1")
    expect_identical(captured_call$geom.text.size, 6)
    expect_identical(colnames(captured_call$data), c("S1", "S2", "S3", "S4"))
    expect_identical(rownames(captured_call$data), assay_data$Name)
    expect_identical(captured_call$design_matrix$Run, design_matrix$Run)
})

test_that("plotRle resolves through the active lipid S4 plotting source", {
    wrapper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_objects.R"))
    helper_path <- file.path("..", "..", "R", "func_lipid_s4_plotting_methods.R")
    helper_source_lines <- if (file.exists(helper_path)) readLines(helper_path) else character()
    helper_defines_plot_rle <- any(grepl("f = \"plotRle\"", helper_source_lines, fixed = TRUE))

    if (helper_defines_plot_rle) {
        expect_false(any(grepl("f = \"plotRle\"", wrapper_source_lines, fixed = TRUE)))
        expect_true(helper_defines_plot_rle)
    } else {
        expect_true(any(grepl("f = \"plotRle\"", wrapper_source_lines, fixed = TRUE)))
    }
})

test_that("plotRle passes labeled assay matrices and grouping info to the helper", {
    assay_data <- tibble::tibble(
        Name = paste0("L", seq_len(5)),
        Annotation = paste0("A", seq_len(5)),
        S1 = c(10, 11, 12, 13, 14),
        S2 = c(12, 11, 13, 14, 15),
        S3 = c(20, 19, 18, 17, 16),
        S4 = c(18, 17, 16, 15, 14)
    )
    design_matrix <- data.frame(
        Run = c("S1", "S2", "S3", "S4"),
        group = c("A", "A", "B", "B"),
        display_name = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
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

    had_plot_helper <- exists("plotRleHelper", envir = .GlobalEnv, inherits = FALSE)
    old_plot_helper <- if (had_plot_helper) {
        get("plotRleHelper", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        if (had_plot_helper) {
            assign("plotRleHelper", old_plot_helper, envir = .GlobalEnv)
        } else if (exists("plotRleHelper", envir = .GlobalEnv, inherits = FALSE)) {
            rm("plotRleHelper", envir = .GlobalEnv)
        }
    }, add = TRUE)

    captured_call <- NULL
    assign(
        "plotRleHelper",
        function(Y, rowinfo = NULL, yaxis_limit = c(-0.5, 0.5), ...) {
            captured_call <<- list(
                Y = Y,
                rowinfo = rowinfo,
                yaxis_limit = yaxis_limit
            )

            ggplot2::ggplot(
                data.frame(
                    sample = rownames(Y),
                    value = seq_len(nrow(Y))
                ),
                ggplot2::aes(sample, value)
            )
        },
        envir = .GlobalEnv
    )

    rle_plots <- plotRle(
        lipid_object,
        grouping_variable = "group",
        yaxis_limit = c(-1, 1),
        sample_label = "display_name"
    )

    expect_named(rle_plots, "Assay1")
    expect_s3_class(rle_plots$Assay1, "ggplot")
    expect_identical(unname(rownames(captured_call$Y)), design_matrix$display_name)
    expect_identical(colnames(captured_call$Y), assay_data$Name)
    expect_identical(unname(captured_call$rowinfo), design_matrix$group)
    expect_identical(captured_call$yaxis_limit, c(-1, 1))
})

test_that("plotDensity resolves through the active lipid S4 plotting source", {
    wrapper_source_lines <- readLines(file.path("..", "..", "R", "func_lipid_s4_objects.R"))
    helper_path <- file.path("..", "..", "R", "func_lipid_s4_plotting_methods.R")
    helper_source_lines <- if (file.exists(helper_path)) readLines(helper_path) else character()
    helper_plot_density_defs <- sum(grepl("f = \"plotDensity\"", helper_source_lines, fixed = TRUE))

    if (helper_plot_density_defs > 0) {
        expect_false(any(grepl("f = \"plotDensity\"", wrapper_source_lines, fixed = TRUE)))
        expect_gte(helper_plot_density_defs, 2)
    } else {
        expect_identical(sum(grepl("f = \"plotDensity\"", wrapper_source_lines, fixed = TRUE)), 2L)
    }
})

test_that("plotDensity on LipidomicsAssayData returns an empty list for empty assay inputs", {
    empty_object <- createLipidomicsAssayData(
        lipid_data = list(),
        design_matrix = data.frame(
            Run = "S1",
            group = "A",
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "Name",
        annotation_id_column = "Annotation",
        sample_id = "Run",
        group_id = "group"
    )

    expect_warning(
        density_plots <- plotDensity(
            empty_object,
            grouping_variable = "group",
            title = "QC density"
        ),
        "No assays found in `lipid_data` slot"
    )

    expect_identical(density_plots, list())
})

test_that("plotDensity converts PCA ggplot input into titled density patchworks", {
    pca_plots <- list(
        Assay1 = ggplot2::ggplot(
            tibble::tibble(
                PC1 = c(1, 1.5, -1, -1.5),
                PC2 = c(2, 2.5, -2, -2.5),
                group = c("A", "A", "B", "B")
            ),
            ggplot2::aes(PC1, PC2, colour = group)
        ) +
            ggplot2::labs(title = "Original PCA")
    )

    density_plots <- plotDensity(
        pca_plots,
        grouping_variable = "group",
        title = "QC density",
        font_size = 6
    )

    expect_named(density_plots, "Assay1")
    expect_true(any(inherits(density_plots$Assay1, c("patchwork", "ggplot"))))
    expect_identical(
        density_plots$Assay1$patches$plots[[1]]$labels$title,
        "QC density - Assay1"
    )
    expect_true(all(c("PC1", "PC2", "group") %in% colnames(density_plots$Assay1$patches$plots[[1]]$data)))
    expect_equal(
        sort(unique(as.character(density_plots$Assay1$patches$plots[[1]]$data$group))),
        c("A", "B")
    )
    expect_identical(
        density_plots$Assay1$patches$plots[[1]]$theme$text$size,
        6
    )
})
