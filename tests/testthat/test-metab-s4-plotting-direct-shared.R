# fidelity-coverage-compare: shared
library(testthat)
library(ggplot2)

createMetaboliteAssayData <- get(
  "createMetaboliteAssayData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

localNamespaceBinding <- function(name, value, package = "MultiScholaR", env = parent.frame()) {
  ns <- asNamespace(package)
  had_value <- exists(name, envir = ns, inherits = FALSE)
  old_value <- if (had_value) get(name, envir = ns, inherits = FALSE) else NULL
  was_locked <- had_value && bindingIsLocked(name, ns)

  if (was_locked) {
    unlockBinding(name, ns)
  }
  assign(name, value, envir = ns)
  if (was_locked) {
    lockBinding(name, ns)
  }

  withr::defer({
    if (exists(name, envir = ns, inherits = FALSE) && bindingIsLocked(name, ns)) {
      unlockBinding(name, ns)
    }
    if (had_value) {
      assign(name, old_value, envir = ns)
    } else if (exists(name, envir = ns, inherits = FALSE)) {
      rm(list = name, envir = ns)
    }
    if (was_locked) {
      lockBinding(name, ns)
    }
  }, envir = env)
}

newDirectMetabPlotObject <- function(named = TRUE, sample_ids = c("S1", "S2", "S3")) {
  assay_tbl <- data.frame(
    Name = c("M1", "M2", "M3", "M4"),
    annotation = c("alpha", "beta", "gamma", "delta"),
    S1 = c(10, 22, 35, 47),
    S2 = c(12, 19, 32, 49),
    S3 = c(40, 51, 67, 82),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  assay_list <- list(assay_tbl)
  if (named) {
    names(assay_list) <- "LCMS_Pos"
  }

  createMetaboliteAssayData(
    metabolite_data = assay_list,
    design_matrix = data.frame(
      Run = sample_ids,
      Group = c("Control", "Control", "Treatment"),
      Batch = c("B1", "B1", "B2"),
      Label = c("Sample 1", "Sample 2", "Sample 3"),
      stringsAsFactors = FALSE
    ),
    metabolite_id_column = "Name",
    annotation_id_column = "annotation",
    sample_id = "Run",
    group_id = "Group"
  )
}

test_that("metabolomics S4 plotPca preserves helper delegation, title composition, and unnamed assay fallback", {
  captured <- new.env(parent = emptyenv())
  localNamespaceBinding(
    "plotPcaHelper",
    function(data,
             design_matrix,
             sample_id_column,
             grouping_variable,
             shape_variable = NULL,
             label_column = NULL,
             title,
             geom.text.size = 8,
             ...) {
      captured$call <- list(
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
          Group = c("Control", "Treatment")
        ),
        ggplot2::aes(PC1, PC2, colour = Group)
      ) +
        ggplot2::labs(title = title)
    }
  )

  plots <- plotPca(
    newDirectMetabPlotObject(),
    grouping_variable = "Group",
    shape_variable = "Batch",
    label_column = "Label",
    title = "QC PCA",
    font_size = 6
  )

  expect_named(plots, "LCMS_Pos")
  expect_s3_class(plots$LCMS_Pos, "ggplot")
  expect_identical(plots$LCMS_Pos$labels$title, "QC PCA - LCMS_Pos")
  expect_identical(captured$call$sample_id_column, "Run")
  expect_identical(captured$call$grouping_variable, "Group")
  expect_identical(captured$call$shape_variable, "Batch")
  expect_identical(captured$call$label_column, "Label")
  expect_identical(captured$call$title, "QC PCA - LCMS_Pos")
  expect_identical(captured$call$geom.text.size, 6)
  expect_identical(rownames(captured$call$data), c("M1", "M2", "M3", "M4"))
  expect_identical(colnames(captured$call$data), c("S1", "S2", "S3"))

  expect_warning(
    unnamed_plots <- plotPca(
      newDirectMetabPlotObject(named = FALSE),
      grouping_variable = "Group"
    ),
    "Assay list was unnamed. Using default names",
    fixed = TRUE
  )
  expect_named(unnamed_plots, "Assay_1")

  expect_error(
    plotPca(newDirectMetabPlotObject(), grouping_variable = "MissingGroup"),
    "`grouping_variable` 'MissingGroup' not found in design_matrix.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 plotRle preserves helper delegation, label remapping, and validation", {
  captured <- new.env(parent = emptyenv())
  localNamespaceBinding(
    "plotRleHelper",
    function(Y, rowinfo = NULL, yaxis_limit = c(), ...) {
      captured$call <- list(
        Y = Y,
        rowinfo = rowinfo,
        yaxis_limit = yaxis_limit
      )

      ggplot2::ggplot(
        data.frame(sample = rownames(Y), value = seq_len(nrow(Y))),
        ggplot2::aes(sample, value)
      )
    }
  )

  expect_message(
    rle_plots <- plotRle(
      newDirectMetabPlotObject(),
      grouping_variable = "Group",
      yaxis_limit = c(-2, 2),
      sample_label = "Label"
    ),
    "--- Entering plotRle for MetaboliteAssayData ---"
  )

  expect_named(rle_plots, "LCMS_Pos")
  expect_s3_class(rle_plots$LCMS_Pos, "ggplot")
  expect_identical(unname(rownames(captured$call$Y)), c("Sample 1", "Sample 2", "Sample 3"))
  expect_identical(colnames(captured$call$Y), c("M1", "M2", "M3", "M4"))
  expect_identical(unname(captured$call$rowinfo), c("Control", "Control", "Treatment"))
  expect_identical(captured$call$yaxis_limit, c(-2, 2))

  expect_warning(
    unnamed_rle <- plotRle(
      newDirectMetabPlotObject(named = FALSE),
      grouping_variable = "Group"
    ),
    "Assay list was unnamed. Using default names",
    fixed = TRUE
  )
  expect_named(unnamed_rle, "Assay_1")

  expect_error(
    plotRle(
      newDirectMetabPlotObject(),
      grouping_variable = "Group",
      sample_label = "MissingLabel"
    ),
    "`sample_label` 'MissingLabel' not found in design_matrix.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 plotDensity preserves list handling and object PCA branches", {
  pca_plot <- ggplot2::ggplot(
    data.frame(
      PC1 = c(1, 2, -1),
      PC2 = c(2, 3, -2),
      Group = c("A", "A", "B")
    ),
    ggplot2::aes(PC1, PC2, colour = Group)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Original PCA")

  list_plots <- plotDensity(
    list(LCMS_Pos = pca_plot),
    grouping_variable = "Group",
    title = "QC Density",
    font_size = 11
  )
  expect_named(list_plots, "LCMS_Pos")
  expect_true(any(inherits(list_plots$LCMS_Pos, c("patchwork", "ggplot"))))
  expect_identical(list_plots$LCMS_Pos[[1]]$labels$title, "QC Density - LCMS_Pos")
  expect_identical(list_plots$LCMS_Pos[[1]]$labels$y, "PC1")
  expect_identical(list_plots$LCMS_Pos[[2]]$labels$y, "PC2")

  expect_warning(
    unnamed_list_plots <- plotDensity(
      list(pca_plot),
      grouping_variable = "Group",
      title = "QC Density"
    ),
    "Input ggplot list was unnamed. Using default names",
    fixed = TRUE
  )
  expect_named(unnamed_list_plots, "Plot_1")

  if (requireNamespace("mixOmics", quietly = TRUE)) {
    captured <- new.env(parent = emptyenv())
    localNamespaceBinding(
      "pca",
      function(X, ncomp = 2, ...) {
        captured$X <- X
        list(
          variates = list(
            X = structure(
              matrix(
                c(1, 2, -1, -2, 0.5, 1.5),
                ncol = 2,
                dimnames = list(rownames(X), c("PC1", "PC2"))
              ),
              class = "matrix"
            )
          )
        )
      },
      package = "mixOmics"
    )

    object_density <- plotDensity(
      newDirectMetabPlotObject(),
      grouping_variable = "Group",
      title = "QC Density",
      font_size = 7
    )

    expect_named(object_density, "LCMS_Pos")
    expect_true(any(inherits(object_density$LCMS_Pos, c("patchwork", "ggplot"))))
    expect_identical(colnames(captured$X), c("M1", "M2", "M3", "M4"))
  } else {
    expect_warning(
      expect_warning(
        object_density <- plotDensity(
          newDirectMetabPlotObject(),
          grouping_variable = "Group",
          title = "QC Density",
          font_size = 7
        ),
        "Error during PCA calculation for density plot",
        fixed = TRUE
      ),
        "PCA result is invalid or missing PC1/PC2. Skipping Density plot.",
        fixed = TRUE
      )
    expect_length(object_density, 0)
  }

  expect_warning(
    empty_density <- plotDensity(
      createMetaboliteAssayData(
        metabolite_data = list(),
        design_matrix = data.frame(
          Run = "S1",
          Group = "A",
          stringsAsFactors = FALSE
        ),
        metabolite_id_column = "Name",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "Group"
      ),
      grouping_variable = "Group",
      title = "Empty Density"
    ),
    "No assays found in `metabolite_data` slot",
    fixed = TRUE
  )
  expect_identical(empty_density, list())
})
