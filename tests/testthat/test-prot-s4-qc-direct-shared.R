# fidelity-coverage-compare: shared
library(testthat)

msr <- function(name) {
  get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

ProteinQuantitativeData <- msr("ProteinQuantitativeData")
removeProteinsWithOnlyOneReplicate <- msr("removeProteinsWithOnlyOneReplicate")
removeRowsWithMissingValuesPercent <- msr("removeRowsWithMissingValuesPercent")
filterSamplesByProteinCorrelationThreshold <- msr("filterSamplesByProteinCorrelationThreshold")
cleanDesignMatrix <- msr("cleanDesignMatrix")
proteinIntensityFiltering <- msr("proteinIntensityFiltering")
plotRle <- msr("plotRle")
plotRleList <- msr("plotRleList")
savePlotRleList <- msr("savePlotRleList")
plotPca <- msr("plotPca")
plotPcaList <- msr("plotPcaList")
savePlotPcaList <- msr("savePlotPcaList")
getPcaMatrix <- msr("getPcaMatrix")
proteinTechRepCorrelation <- msr("proteinTechRepCorrelation")
plotPearson <- msr("plotPearson")
pearsonCorForSamplePairs <- msr("pearsonCorForSamplePairs")
plotPcaBox <- msr("plotPcaBox")
plotDensityList <- msr("plotDensityList")
savePlotDensityList <- msr("savePlotDensityList")
filterMinNumPeptidesPerProtein <- msr("filterMinNumPeptidesPerProtein")

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

newProtS4QcObject <- function(args = list(), protein_quant_table = NULL, design_matrix = NULL) {
  if (is.null(protein_quant_table)) {
    protein_quant_table <- data.frame(
      Protein.Ids = c("P1", "P2", "P3"),
      S1 = c(10, 20, 30),
      S2 = c(11, 19, 31),
      S3 = c(5, 6, 7),
      S4 = c(6, 7, 8),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("Control", "Control", "Pool", "Pool"),
      batch = c("B1", "B2", "B1", "B2"),
      replicates = c("rep1", "rep1", "pool_pair", "pool_pair"),
      label = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
      stringsAsFactors = FALSE
    )
  }

  ProteinQuantitativeData(
    protein_quant_table = protein_quant_table,
    design_matrix = design_matrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    protein_id_column = "Protein.Ids",
    args = args
  )
}

test_that("protein S4 helper-entry methods preserve current registration and failure behavior", {
  object <- newProtS4QcObject()

  if (methods::hasMethod("removeProteinsWithOnlyOneReplicate", "ProteinQuantitativeData")) {
    expect_error(
      removeProteinsWithOnlyOneReplicate(
        object,
        core_utilisation = 1,
        grouping_variable = "group"
      ),
      "invalid argument type"
    )
  } else {
    expect_false(methods::hasMethod("removeProteinsWithOnlyOneReplicate", "ProteinQuantitativeData"))
  }

  if (methods::hasMethod("removeRowsWithMissingValuesPercent", "ProteinQuantitativeData")) {
    expect_error(
      suppressMessages(
        removeRowsWithMissingValuesPercent(
          object,
          ruv_grouping_variable = "group",
          groupwise_percentage_cutoff = 50,
          max_groups_percentage_cutoff = 50,
          proteins_intensity_cutoff_percentile = 0.5
        )
      ),
      "invalid argument type"
    )
  } else {
    expect_false(methods::hasMethod("removeRowsWithMissingValuesPercent", "ProteinQuantitativeData"))
  }

  if (methods::hasMethod("filterSamplesByProteinCorrelationThreshold", "ProteinQuantitativeData")) {
    expect_error(
      suppressMessages(
        filterSamplesByProteinCorrelationThreshold(
          object,
          pearson_correlation_per_pair = data.frame(
            Run.x = "S1",
            Run.y = "S2",
            pearson_correlation = 0.99,
            stringsAsFactors = FALSE
          ),
          min_pearson_correlation_threshold = 0.75
        )
      ),
      "invalid argument type"
    )
  } else {
    expect_false(methods::hasMethod("filterSamplesByProteinCorrelationThreshold", "ProteinQuantitativeData"))
  }

  if (methods::hasMethod("proteinIntensityFiltering", "ProteinQuantitativeData")) {
    expect_error(
      proteinIntensityFiltering(
        object,
        proteins_intensity_cutoff_percentile = 0.5,
        proteins_proportion_of_samples_below_cutoff = 0.5,
        core_utilisation = 1
      ),
      "first argument must be a vector"
    )
  } else {
    expect_false(methods::hasMethod("proteinIntensityFiltering", "ProteinQuantitativeData"))
  }
})

test_that("cleanDesignMatrix preserves sample ordering and current warning fallbacks", {
  reordered_object <- newProtS4QcObject(
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2"),
      S2 = c(20, 21),
      S1 = c(10, 11),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("Control", "Treatment"),
      stringsAsFactors = FALSE
    )
  )

  cleaned <- cleanDesignMatrix(reordered_object)
  expect_identical(cleaned@design_matrix$Run, c("S2", "S1"))
  expect_identical(cleaned@design_matrix$group, c("Treatment", "Control"))

  missing_sample_column <- reordered_object
  missing_sample_column@sample_id <- "MissingRun"
  expect_warning(
    unchanged_missing <- cleanDesignMatrix(missing_sample_column),
    "Could not extract valid sample IDs"
  )
  expect_identical(unchanged_missing@design_matrix$Run, c("S1", "S2"))

  unmatched_object <- reordered_object
  colnames(unmatched_object@protein_quant_table) <- c("Protein.Ids", "Z1", "Z2")
  expect_warning(
    unchanged_unmatched <- cleanDesignMatrix(unmatched_object),
    "No sample columns identified"
  )
  expect_identical(unchanged_unmatched@design_matrix$Run, c("S1", "S2"))
})

test_that("plotRle methods preserve sample relabeling, row grouping, and file naming", {
  object <- newProtS4QcObject()
  object@protein_quant_table$S1[[1]] <- Inf
  helper_calls <- list()
  saved_files <- character()

  testthat::local_mocked_bindings(
    plotRleHelper = function(input_matrix, rowinfo, yaxis_limit) {
      helper_calls[[length(helper_calls) + 1L]] <<- list(
        input_matrix = input_matrix,
        rowinfo = rowinfo,
        yaxis_limit = yaxis_limit
      )
      ggplot2::ggplot(
        data.frame(
          x = seq_along(rowinfo),
          y = seq_along(rowinfo),
          rowinfo = rowinfo,
          stringsAsFactors = FALSE
        ),
        ggplot2::aes(x, y, color = rowinfo)
      ) +
        ggplot2::geom_point()
    },
    ggsave = function(plot, filename) {
      saved_files <<- c(saved_files, filename)
      invisible(filename)
    },
    .package = "MultiScholaR"
  )

  rle_plot <- plotRle(
    object,
    grouping_variable = "group",
    yaxis_limit = c(-2, 2),
    sample_label = "label"
  )
  rle_list <- plotRleList(object, list_of_columns = c("group", "batch"), yaxis_limit = c(-1, 1))

  expect_s3_class(rle_plot, "ggplot")
  expect_equal(helper_calls[[1L]]$rowinfo, c("Control", "Control", "Pool", "Pool"))
  expect_true(anyNA(helper_calls[[1L]]$input_matrix[, 1]))
  expect_named(rle_list, c("group", "batch"))
  expect_equal(helper_calls[[2L]]$rowinfo, c("Control", "Control", "Pool", "Pool"))
  expect_equal(helper_calls[[3L]]$rowinfo, c("B1", "B2", "B1", "B2"))

  saved_manifest <- savePlotRleList(
    list(group = rle_plot, batch = rle_list$batch),
    suffix = c("png", "pdf"),
    output_dir = tempdir()
  )

  expect_equal(
    saved_manifest$filename,
    c("RLE_group.png", "RLE_group.pdf", "RLE_batch.png", "RLE_batch.pdf")
  )
  expect_true(all(file.path(tempdir(), saved_manifest$filename) %in% saved_files))
})

test_that("plotPca methods preserve validation, helper payloads, and current file naming", {
  object <- newProtS4QcObject()
  helper_calls <- list()
  saved_files <- character()

  testthat::local_mocked_bindings(
    plotPcaHelper = function(frozen_protein_matrix_pca,
                             design_matrix,
                             sample_id_column,
                             grouping_variable,
                             shape_variable = NULL,
                             label_column = "",
                             title = "",
                             geom.text.size = 8) {
      helper_calls[[length(helper_calls) + 1L]] <<- list(
        matrix = frozen_protein_matrix_pca,
        design_matrix = design_matrix,
        sample_id_column = sample_id_column,
        grouping_variable = grouping_variable,
        shape_variable = shape_variable,
        label_column = label_column,
        title = title,
        geom.text.size = geom.text.size
      )

      plot_data <- data.frame(
        Run = design_matrix[[sample_id_column]],
        PC1 = seq_len(nrow(design_matrix)),
        PC2 = rev(seq_len(nrow(design_matrix))),
        stringsAsFactors = FALSE
      )
      plot_data[[grouping_variable]] <- design_matrix[[grouping_variable]]
      if (!is.null(shape_variable)) {
        plot_data[[shape_variable]] <- design_matrix[[shape_variable]]
      }
      if (!identical(label_column, "")) {
        plot_data[[label_column]] <- design_matrix[[label_column]]
      }

      ggplot2::ggplot(plot_data, ggplot2::aes(PC1, PC2, color = .data[[grouping_variable]])) +
        ggplot2::geom_point() +
        ggplot2::labs(title = title)
    },
    plotPcaListHelper = function(frozen_protein_matrix_pca,
                                 design_matrix,
                                 sample_id_column,
                                 grouping_variables_list,
                                 label_column = "",
                                 title = "",
                                 geom.text.size = 8) {
      helper_calls[[length(helper_calls) + 1L]] <<- list(
        grouping_variables_list = grouping_variables_list,
        label_column = label_column,
        title = title,
        geom.text.size = geom.text.size
      )

      stats::setNames(
        lapply(grouping_variables_list, function(group_var) {
          ggplot2::ggplot(
            data.frame(
              PC1 = 1:2,
              PC2 = 2:1,
              group = c(group_var, group_var),
              stringsAsFactors = FALSE
            ),
            ggplot2::aes(PC1, PC2, color = group)
          ) +
            ggplot2::geom_point()
        }),
        grouping_variables_list
      )
    },
    ggsave = function(plot, filename) {
      saved_files <<- c(saved_files, filename)
      invisible(filename)
    },
    .package = "MultiScholaR"
  )

  pca_plot <- plotPca(
    object,
    grouping_variable = "group",
    shape_variable = "batch",
    label_column = "label",
    title = "Protein PCA",
    font_size = 9
  )

  expect_s3_class(pca_plot, "ggplot")
  expect_identical(helper_calls[[1L]]$sample_id_column, "Run")
  expect_identical(helper_calls[[1L]]$grouping_variable, "group")
  expect_identical(helper_calls[[1L]]$shape_variable, "batch")
  expect_identical(helper_calls[[1L]]$label_column, "label")
  expect_identical(helper_calls[[1L]]$title, "Protein PCA")

  expect_error(
    plotPca(object, grouping_variable = 1, label_column = "", title = "bad"),
    "single character string"
  )
  expect_error(
    plotPca(object, grouping_variable = "missing", label_column = "", title = "bad"),
    "not found in design matrix"
  )
  expect_error(
    plotPca(
      object,
      grouping_variable = "group",
      shape_variable = "missing",
      label_column = "",
      title = "bad"
    ),
    "shape_variable 'missing' not found"
  )

  pca_list <- plotPcaList(
    object,
    grouping_variables_list = c("group", "batch"),
    label_column = "label",
    title = "Protein PCA list",
    font_size = 7
  )

  expect_named(pca_list, c("group", "batch"))
  expect_identical(helper_calls[[2L]]$grouping_variables_list, c("group", "batch"))

  saved_manifest <- savePlotPcaList(
    pca_list,
    suffix = c("png", "pdf"),
    output_dir = tempdir()
  )

  expect_equal(
    saved_manifest$filename,
    c("RLE_group.png", "RLE_group.pdf", "RLE_batch.png", "RLE_batch.pdf")
  )
  expect_true(all(file.path(tempdir(), saved_manifest$filename) %in% saved_files))
})

test_that("getPcaMatrix preserves projected sample annotations", {
  skip_if_not_installed("mixOmics")

  object <- newProtS4QcObject(
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2", "P3", "P4"),
      S1 = c(1, 2, 3, 4),
      S2 = c(2, 4, 6, 8),
      S3 = c(8, 6, 4, 2),
      S4 = c(4, 3, 2, 1),
      stringsAsFactors = FALSE
    )
  )

  pca_matrix <- getPcaMatrix(object)

  expect_true(all(c("Run", "group", "batch", "replicates") %in% names(pca_matrix)))
  expect_identical(sort(pca_matrix$Run), c("S1", "S2", "S3", "S4"))
})

test_that("protein correlation and density helpers preserve grouping, filtering, and downstream plotting", {
  object <- newProtS4QcObject()
  captured <- new.env(parent = emptyenv())
  captured$tech_rep_args <- NULL
  captured$density_calls <- list()
  saved_files <- character()

  testthat::local_mocked_bindings(
    proteinTechRepCorrelationHelper = function(design_matrix,
                                               frozen_protein_matrix_pca,
                                               protein_id_column,
                                               sample_id_column,
                                               tech_rep_column,
                                               tech_rep_num_column = NULL,
                                               tech_rep_remove_regex = NULL) {
      captured$tech_rep_args <- list(
        design_matrix = design_matrix,
        matrix = frozen_protein_matrix_pca,
        protein_id_column = protein_id_column,
        sample_id_column = sample_id_column,
        tech_rep_column = tech_rep_column,
        tech_rep_num_column = tech_rep_num_column,
        tech_rep_remove_regex = tech_rep_remove_regex
      )
      data.frame(pair = "rep1", pearson = 0.99, stringsAsFactors = FALSE)
    },
    plotDensity = function(theObject, grouping_variable, title = "", font_size = 8) {
      captured$density_calls[[length(captured$density_calls) + 1L]] <- list(
        grouping_variable = grouping_variable,
        title = title,
        font_size = font_size
      )
      if (identical(grouping_variable, "bad")) {
        stop("bad grouping")
      }
      ggplot2::ggplot(
        data.frame(x = 1:2, y = c(1, 2), group = grouping_variable, stringsAsFactors = FALSE),
        ggplot2::aes(x, y)
      ) +
        ggplot2::geom_line()
    },
    ggsave = function(plot, filename) {
      saved_files <<- c(saved_files, filename)
      invisible(filename)
    },
    .package = "MultiScholaR"
  )

  tech_rep <- proteinTechRepCorrelation(
    object,
    tech_rep_num_column = "pair_index",
    tech_rep_remove_regex = "pool"
  )
  expect_identical(tech_rep$pair, "rep1")
  expect_identical(captured$tech_rep_args$tech_rep_num_column, "pair_index")
  expect_identical(captured$tech_rep_args$tech_rep_remove_regex, "pool")

  all_pairs <- suppressMessages(pearsonCorForSamplePairs(object, exclude_pool_samples = FALSE))
  default_pairs <- suppressMessages(pearsonCorForSamplePairs(object))
  override_pairs <- suppressMessages(pearsonCorForSamplePairs(object, correlation_group = "group"))

  expect_equal(nrow(all_pairs), 2L)
  expect_equal(nrow(default_pairs), 1L)
  expect_identical(default_pairs$Run.x, "S1")
  expect_identical(default_pairs$Run.y, "S2")
  expect_identical(default_pairs$replicates, "rep1")
  expect_equal(nrow(override_pairs), 1L)
  expect_identical(override_pairs$group, "Control")

  no_match_object <- object
  no_match_object@design_matrix$Run <- paste0("X", seq_len(nrow(no_match_object@design_matrix)))
  expect_error(
    suppressMessages(pearsonCorForSamplePairs(no_match_object)),
    "No matching samples found"
  )

  pearson_plot <- suppressMessages(plotPearson(object))
  expect_s3_class(pearson_plot, "ggplot")
  expect_identical(pearson_plot$labels$x, "Pearson Correlation")

  expect_warning(
    density_list <- plotDensityList(
      object,
      grouping_variables_list = c("group", "bad", "batch"),
      title = "Density",
      font_size = 9
    ),
    "Error creating density plot for bad"
  )
  expect_named(density_list, c("group", "batch"))
  expect_equal(vapply(captured$density_calls, `[[`, character(1), "grouping_variable"), c("group", "bad", "batch"))

  saved_manifest <- savePlotDensityList(density_list, suffix = c("png", "pdf"), output_dir = tempdir())
  expect_equal(
    saved_manifest$filename,
    c("Density_group.png", "Density_group.pdf", "Density_batch.png", "Density_batch.pdf")
  )
  expect_true(all(file.path(tempdir(), saved_manifest$filename) %in% saved_files))
})

test_that("plotPcaBox preserves ggplot extraction paths and current validation behavior", {
  skip_if_not_installed("patchwork")

  direct_plot <- ggplot2::ggplot(
    data.frame(
      PC1 = c(1, 2, 3, 4),
      PC2 = c(4, 3, 2, 1),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    ggplot2::aes(PC1, PC2, color = group)
  ) +
    ggplot2::geom_point()

  direct_box <- plotPcaBox(direct_plot, grouping_variable = "group", title = "PCA box", show_legend = TRUE)
  expect_false(is.null(direct_box))

  env_box <- local({
    data <- data.frame(
      PC1 = c(1, 2, 3, 4),
      PC2 = c(4, 3, 2, 1),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
    plot <- ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, color = group)) +
      ggplot2::geom_point()
    plot$data <- NULL
    plotPcaBox(plot, grouping_variable = "group")
  })
  expect_false(is.null(env_box))

  expect_error(plotPcaBox(list(), grouping_variable = "group"), "ggplot object")
  expect_error(
    plotPcaBox(
      ggplot2::ggplot(data.frame(PC1 = 1, PC2 = 2), ggplot2::aes(PC1, PC2)) + ggplot2::geom_point(),
      grouping_variable = "group"
    ),
    "not found in the data"
  )
})

test_that("filterMinNumPeptidesPerProtein preserves peptide evidence filtering and EList synchronization", {
  skip_if_not_installed("limma")

  peptide_summary <- data.frame(
    Protein.Ids = c("P1", "P2", "P3"),
    peptide_count = c(1, 3, 4),
    peptidoform_count = c(1, 2, 5),
    stringsAsFactors = FALSE
  )

  quantified_elist <- structure(
    list(
      E = matrix(1:12, nrow = 3, byrow = TRUE),
      genes = data.frame(protein.id = c("P1", "P2", "P3"), stringsAsFactors = FALSE)
    ),
    class = "EList"
  )

  object <- newProtS4QcObject(
    args = list(
      limpa_dpc_quant_results = list(
        peptide_counts_per_protein = peptide_summary,
        quantified_elist = quantified_elist
      )
    )
  )

  filtered <- filterMinNumPeptidesPerProtein(
    object,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    verbose = TRUE
  )

  expect_identical(filtered@protein_quant_table$Protein.Ids, c("P2", "P3"))
  expect_identical(filtered@args$limpa_dpc_quant_results$quantified_elist$genes$protein.id, c("P2", "P3"))

  missing_summary_object <- newProtS4QcObject(args = list(limpa_dpc_quant_results = list()))
  expect_error(
    filterMinNumPeptidesPerProtein(
      missing_summary_object,
      num_peptides_per_protein_thresh = 2,
      num_peptidoforms_per_protein_thresh = 2,
      verbose = FALSE
    ),
    "Could not find the peptide summary table"
  )
})
