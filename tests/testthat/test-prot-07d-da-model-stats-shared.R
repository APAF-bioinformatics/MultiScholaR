# fidelity-coverage-compare: shared
library(testthat)

localProtDaNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withProtDaNamespaceBindings <- function(replacements, code) {
  mock_frame <- parent.frame()
  package_ns <- asNamespace("MultiScholaR")
  imports_env <- parent.env(package_ns)

  for (name in names(replacements)) {
    target_env <- if (exists(name, envir = package_ns, inherits = FALSE)) package_ns else imports_env
    localProtDaNamespaceBinding(target_env, name, replacements[[name]], mock_frame)
  }

  force(code)
}

is.na.prot_da_rows <- function(x) {
  FALSE
}

is.na.prot_da_weights <- function(x) {
  FALSE
}

names.prot_da_rows <- function(x) {
  c("A", "B")
}

`[[.prot_da_rows` <- function(x, i, ...) {
  attr(x, i, exact = TRUE)
}

`[.prot_da_weights` <- function(x, i, j, ..., drop = FALSE) {
  matrix(1, nrow = length(i), ncol = 1, dimnames = list(i, "weight"))
}

makeProtDaModelObject <- function(args = list()) {
  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      S1 = c(10, 11, 12),
      S2 = c(12, 13, 14),
      S3 = c(20, 21, 22),
      S4 = c(22, 23, 24),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      gene_names = c("GENE1", "GENE2", "GENE3"),
      stringsAsFactors = FALSE
    ),
    args = args
  )
}

makeProtDaMetaboliteObject <- function(args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      assay = data.frame(
        uniprot_acc = c("P1", "P2", "P3"),
        S1 = c(10, 11, 12),
        S2 = c(12, 13, 14),
        S3 = c(20, 21, 22),
        S4 = c(22, 23, 24),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "uniprot_acc",
    annotation_id_column = "uniprot_acc",
    database_identifier_type = "test",
    internal_standard_regex = NA_character_,
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = args
  )
}

makeProtDaMatrix <- function(rows = 4, cols = 4) {
  values <- seq_len(rows * cols)
  matrix(
    values,
    nrow = rows,
    dimnames = list(paste0("P", seq_len(rows)), paste0("S", seq_len(cols)))
  )
}

makeProtDaDesign <- function(include_replicates = FALSE, duplicated_replicates = FALSE) {
  design <- data.frame(
    Run = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  if (include_replicates) {
    design$replicates <- if (duplicated_replicates) c(1, 1, 1, 1) else c(1, 2, 1, 2)
  }
  rownames(design) <- design$Run
  design
}

makeProtDaFitMocks <- function(captured = new.env(parent = emptyenv())) {
  captured$make_contrasts <- list()
  captured$lmfit_calls <- list()
  captured$duplicate_correlation <- 0L

  list(
    makeContrasts = function(contrasts, levels) {
      captured$make_contrasts[[length(captured$make_contrasts) + 1L]] <<- list(
        contrasts = contrasts,
        levels = levels
      )
      matrix(seq_along(levels), nrow = length(levels), dimnames = list(levels, contrasts))
    },
    lmFit = function(object, design, block = NULL, correlation = NULL, ...) {
      captured$lmfit_calls[[length(captured$lmfit_calls) + 1L]] <<- list(
        object = object,
        design = design,
        block = block,
        correlation = correlation
      )
      list(object = object, design = design, block = block, correlation = correlation)
    },
    contrasts.fit = function(fit, contrasts) {
      fit$contrasts <- contrasts
      fit
    },
    eBayes = function(fit, trend = FALSE, robust = FALSE, ...) {
      n <- nrow(fit$object)
      if (is.null(n) || is.na(n)) {
        n <- 4L
      }
      list(
        coefficients = matrix(seq_len(n) / 10, ncol = 1),
        df.residual = rep(6, n),
        df.prior = 4,
        s2.prior = 1.5,
        sigma = rep(2, n),
        s2.post = rep(3, n),
        stdev.unscaled = matrix(rep(1, n), ncol = 1),
        t = matrix(seq_len(n) / 5, ncol = 1),
        p.value = matrix(seq(0.01, 0.04, length.out = n), ncol = 1),
        trend = trend,
        robust = robust
      )
    },
    qvalue = function(p) {
      list(q = pmin(1, as.numeric(p) * 2))
    },
    duplicateCorrelation = function(object, design, block, ...) {
      captured$duplicate_correlation <<- captured$duplicate_correlation + 1L
      list(consensus.correlation = 0.25)
    },
    topTable = function(fit, coef, n = Inf, ...) {
      data.frame(
        logFC = c(1.1, -0.8, 0.3, 0.1),
        AveExpr = seq_len(4),
        t = c(4, -3, 2, 1),
        P.Value = c(0.01, 0.02, 0.03, 0.04),
        adj.P.Val = c(0.02, 0.03, 0.04, 0.05),
        B = seq(0.4, 0.1, length.out = 4),
        row.names = paste0("P", 1:4)
      )
    },
    treat = function(fit, lfc) {
      fit$treat_lfc <- lfc
      fit
    },
    topTreat = function(fit, coef, n = Inf, ...) {
      data.frame(
        logFC = c(1.4, -1.2),
        AveExpr = c(3, 4),
        t = c(5, -4),
        P.Value = c(0.001, 0.002),
        adj.P.Val = c(0.002, 0.002),
        B = c(1, 0.8),
        row.names = c("P1", "P2")
      )
    }
  )
}

test_that("proteomics DA model stats helpers execute limma-style analysis paths", {
    skip("quarantined Q-11 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
  captured <- new.env(parent = emptyenv())
  fit_mocks <- makeProtDaFitMocks(captured)
  eb_fit <- getFromNamespace("ebFit", "MultiScholaR")
  run_test <- getFromNamespace("runTest", "MultiScholaR")
  run_tests <- getFromNamespace("runTests", "MultiScholaR")
  run_tests_contrasts <- getFromNamespace("runTestsContrasts", "MultiScholaR")

  withProtDaNamespaceBindings(fit_mocks, {
    data <- makeProtDaMatrix()
    design <- makeProtDaDesign()

    eb_result <- eb_fit(data, design, contr.matrix = matrix(1, nrow = 2))
    expect_s3_class(eb_result$table, "data.frame")
    expect_true(all(c("logFC", "raw_pvalue", "fdr_qvalue", "s2.post") %in% names(eb_result$table)))
    expect_equal(nrow(eb_result$table), nrow(data))

    run_test_result <- run_test(
      ID = rownames(data),
      A = data[, 1:2],
      B = data[, 3:4],
      group_A = "A",
      group_B = "B",
      design_matrix = design,
      formula_string = "~ 0 + group"
    )
    expect_s3_class(run_test_result$table, "data.frame")
    expect_equal(unique(run_test_result$table$comparison), "log(B) minus log(A)")
    expect_equal(run_test_result$table$meanA, unname(rowMeans(data[, 1:2])))
    expect_equal(run_test_result$table$meanB, unname(rowMeans(data[, 3:4])))

    run_tests_result <- run_tests(
      ID = rownames(data),
      data = data,
      test_pairs = data.frame(A = "A", B = "B", stringsAsFactors = FALSE),
      sample_columns = colnames(data),
      sample_rows_list = NA,
      type_of_grouping = list(A = c("S1", "S2"), B = c("S3", "S4")),
      design_matrix = design,
      formula_string = "~ 0 + group"
    )
    expect_named(run_tests_result, "B vs A")
    expect_s3_class(run_tests_result[["B vs A"]]$results, "data.frame")
    expect_equal(nrow(run_tests_result[["B vs A"]]$counts), 4)

    row_filter <- list(A = c("P1", "P2"))
    withProtDaNamespaceBindings(
      list(
        runTest = function(ID, A, B, group_A, group_B, design_matrix, formula_string,
                           contrast_variable = "group", weights = NA) {
          captured$filtered_rows <<- rownames(A)
          list(table = data.frame(marker = 1, stringsAsFactors = FALSE), fit.eb = "fit")
        }
      ),
      {
        filtered_result <- run_tests(
          ID = rownames(data),
          data = data,
          test_pairs = data.frame(A = "A", B = "B", stringsAsFactors = FALSE),
          sample_columns = colnames(data),
          sample_rows_list = row_filter,
          type_of_grouping = list(A = c("S1", "S2"), B = c("S3", "S4")),
          design_matrix = design,
          formula_string = "~ 0 + group"
        )
      }
    )
    expect_named(filtered_result, "B vs A")
    expect_equal(captured$filtered_rows, c("P1", "P2"))

    contrast_result <- run_tests_contrasts(
      data = data,
      contrast_strings = "A_vs_B=groupB-groupA",
      design_matrix = design,
      formula_string = "~ 0 + group",
      treat_lfc_cutoff = NA
    )
    expect_named(contrast_result$results, "A_vs_B=groupB-groupA")
    expect_true(all(c("raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment") %in% names(contrast_result$results[[1]])))
    expect_length(contrast_result$qvalue_warnings, 0)

    # KIT-064 addendum (roar BUG-21): on the edge path (fewer than three p-values, so qvalue() is not attempted) the
    # contrast's q-value warning is recorded and nothing leaks into the global environment; the former `<<-` in the
    # contrast function's own frame did neither.
    if (exists("qvalue_failed", envir = globalenv(), inherits = FALSE)) rm("qvalue_failed", envir = globalenv())
    edge_result <- run_tests_contrasts(
      data = data[1:2, , drop = FALSE],
      contrast_strings = "A_vs_B=groupB-groupA",
      design_matrix = design,
      formula_string = "~ 0 + group",
      treat_lfc_cutoff = NA
    )
    expect_true(isTRUE(edge_result$qvalue_warnings[["A_vs_B=groupB-groupA"]]))
    expect_false(exists("qvalue_failed", envir = globalenv(), inherits = FALSE))

    replicate_result <- run_tests_contrasts(
      data = data,
      contrast_strings = "A_vs_B=groupB-groupA",
      design_matrix = makeProtDaDesign(include_replicates = TRUE, duplicated_replicates = TRUE),
      formula_string = "~ 0 + group",
      treat_lfc_cutoff = 0.5
    )
    expect_named(replicate_result$results, "A_vs_B=groupB-groupA")
    # KIT-035: a replicate block equal to the group is spanned by the design; standard lmFit, no duplicateCorrelation
    # (limma would have set the intrablock correlation to zero and warned).
    expect_equal(captured$duplicate_correlation, 0L)
    expect_equal(replicate_result$results[[1]]$fdr_qvalue, p.adjust(c(0.001, 0.002), method = "BH"))
  })
})

test_that("daAnalysisWrapperFunction assembles standard DA result artifacts", {
  captured <- new.env(parent = emptyenv())
  captured$contrast_strings <- NULL
  captured$omit_raw_pvalue <- FALSE
  da_wrapper <- getFromNamespace("daAnalysisWrapperFunction", "MultiScholaR")

  params <- list(
    contrasts_tbl = data.frame(full_format = "A_vs_B=groupB-groupA", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group",
    group_id = "group",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0,
    eBayes_trend = FALSE,
    eBayes_robust = FALSE,
    args_group_pattern = "S[0-9]+",
    args_row_id = "uniprot_acc"
  )

  empty_plot <- function() {
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point()
  }

  withProtDaNamespaceBindings(
    list(
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        if (param_name_string %in% names(params)) {
          return(params[[param_name_string]])
        }
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      },
      plotRle = function(...) empty_plot(),
      plotPca = function(...) empty_plot(),
      plotNumOfValuesNoLog = function(...) empty_plot(),
      getDataMatrix = function(theObject) {
        if (inherits(theObject, "MetaboliteAssayData")) {
          return(as.matrix(theObject@metabolite_data[[1]][, c("S1", "S2", "S3", "S4")]))
        }
        as.matrix(theObject@protein_quant_table[, c("S1", "S2", "S3", "S4")])
      },
      runTestsContrasts = function(data, contrast_strings, ...) {
        captured$contrast_strings <<- contrast_strings
        result_table <- data.frame(
          uniprot_acc = c("P1", "P2", "P3"),
          logFC = c(1.2, -0.7, 0.1),
          raw_pvalue = c(0.001, 0.02, 0.5),
          fdr_qvalue = c(0.01, 0.04, 0.8),
          fdr_value_bh_adjustment = c(0.01, 0.04, 0.8),
          stringsAsFactors = FALSE
        )
        if (isTRUE(captured$omit_raw_pvalue)) {
          result_table$raw_pvalue <- NULL
        }
        result_list <- list(result_table)
        names(result_list) <- contrast_strings[[1]]
        list(
          results = result_list,
          fit.eb = "fit"
        )
      },
      getSignificantData = function(...) {
        data.frame(
          uniprot_acc = c("P1", "P2", "P3"),
          comparison = c("A_vs_B", "A_vs_B", "A_vs_B"),
          log_intensity = c("groupB-groupA", "groupB-groupA", "groupB-groupA"),
          raw_pvalue = c(0.001, 0.02, 0.5),
          fdr_qvalue = c(0.01, 0.04, 0.8),
          fdr_value_bh_adjustment = c(0.01, 0.04, 0.8),
          lqm = c(2, 1.4, 0.1),
          logFC = c(1.2, -0.7, 0.1),
          analysis_type = c("RUV applied", "RUV applied", "RUV applied"),
          colour = c("purple", "purple", "black"),
          stringsAsFactors = FALSE
        )
      },
      plotVolcano = function(...) empty_plot(),
      printCountDaGenesTable = function(...) {
        list(
          plot = empty_plot(),
          table = data.frame(
            status = c("Significant and Up", "Significant and Down", "Not significant"),
            counts = c(1, 1, 1),
            comparison = "A_vs_B",
            stringsAsFactors = FALSE
          )
        )
      },
      printPValuesDistribution = function(...) empty_plot(),
      uniprot_tbl = data.frame(
        Entry = c("P1", "P2", "P3"),
        gene_names = c("GENE1; alternate", "", NA_character_),
        stringsAsFactors = FALSE
      ),
      protein_id_column = "Entry",
      printOneVolcanoPlotWithProteinLabel = function(...) empty_plot(),
      getCountsTable = function(theObject) {
        if (inherits(theObject, "MetaboliteAssayData")) {
          return(theObject@metabolite_data[[1]])
        }
        theObject@protein_quant_table
      },
      createDaResultsLongFormat = function(...) {
        data.frame(
          uniprot_acc = c("P1", "P2", "P3"),
          comparison = c("A_vs_B", "A_vs_B", "A_vs_B"),
          raw_pvalue = c(0.001, 0.02, 0.5),
          fdr_qvalue = c(0.01, 0.04, 0.8),
          log2FC = c(1.2, -0.7, 0.1),
          stringsAsFactors = FALSE
        )
      },
      plotOneVolcanoNoVerticalLines = function(...) empty_plot()
    ),
    {
      result <- da_wrapper(makeProtDaModelObject())

      expect_type(result, "list")
      expect_equal(captured$contrast_strings, "A_vs_B=groupB-groupA")
      expect_s3_class(result$rle_plot, "ggplot")
      expect_s3_class(result$pca_plot, "ggplot")
      expect_s3_class(result$da_proteins_wide, "data.frame")
      expect_s3_class(result$da_proteins_long, "data.frame")
      expect_s3_class(result$list_of_volcano_plots_with_gene_names$plot[[1]], "ggplot")
      expect_equal(result$num_sig_da_molecules$counts, c(1L, 1L, 1L))

      params$contrasts_tbl <- data.frame(contrast = "groupB-groupA", stringsAsFactors = FALSE)
      captured$contrast_strings <- NULL
      captured$omit_raw_pvalue <- TRUE
      auto_result <- da_wrapper(makeProtDaModelObject())
      expect_type(auto_result, "list")
      expect_equal(unname(captured$contrast_strings), "B_vs_A=groupB-groupA")
      expect_length(auto_result$raw_pval_histograms, 0)

      captured$omit_raw_pvalue <- FALSE
      expect_error(
        da_wrapper(makeProtDaMetaboliteObject()),
        "`select\\(\\)` doesn't handle lists"
      )
    }
  )
})

test_that("proteomics DA output helpers write expected result artifacts", {
  captured <- new.env(parent = emptyenv())
  captured$vroom_paths <- character()
  captured$xlsx_paths <- character()
  captured$plot_names <- character()
  captured$interactive <- NULL
  output_da_analysis_results <- getFromNamespace("outputDaAnalysisResults", "MultiScholaR")
  save_da_protein_list <- getFromNamespace("saveDaProteinList", "MultiScholaR")

  output_dir <- file.path(tempdir(), "prot-da-output-helper")
  graphs_dir <- file.path(tempdir(), "prot-da-output-graphs")
  unlink(c(output_dir, graphs_dir), recursive = TRUE, force = TRUE)
  dir.create(output_dir, recursive = TRUE)

  empty_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()

  obj <- makeProtDaModelObject(args = list(
    uniprot_tbl = data.frame(
      Entry = c("P1", "P2", "P3"),
      gene_names = c("GENE1 full", "GENE2 full", ""),
      stringsAsFactors = FALSE
    ),
    da_output_dir = output_dir,
    publication_graphs_dir = graphs_dir,
    file_prefix = "da_proteins",
    plots_format = c("pdf"),
    args_row_id = "uniprot_acc",
    da_q_val_thresh = 0.05,
    gene_names_column = "gene_names",
    fdr_column = "fdr_qvalue",
    raw_p_value_column = "raw_pvalue",
    log2fc_column = "log2FC",
    uniprot_id_column = "Entry",
    display_columns = c("gene_name")
  ))

  da_rows <- data.frame(
    uniprot_acc = c("P1-1", "P2-2", "P3-3"),
    comparison = c("A_vs_B", "A_vs_B", "A_vs_B"),
    raw_pvalue = c(0.001, 0.02, 0.5),
    fdr_qvalue = c(0.01, 0.04, 0.8),
    lqm = c(2, 1.4, 0.1),
    log2FC = c(1.2, -0.7, 0.1),
    colour = c("purple", "purple", "black"),
    stringsAsFactors = FALSE
  )

  da_analysis_results <- list(
    theObject = obj,
    pca_plot = empty_plot,
    pca_plot_with_labels = empty_plot,
    rle_plot = empty_plot,
    plot_num_of_values = empty_plot,
    contrasts_results = list(fit.eb = list(sigma = 1)),
    significant_rows = da_rows,
    volplot_plot = empty_plot,
    num_sig_da_molecules_first_go = list(table = data.frame(status = "Up", counts = 1)),
    pvalhist = empty_plot,
    da_proteins_wide = da_rows[, c("uniprot_acc", "comparison", "raw_pvalue", "fdr_qvalue", "log2FC")],
    da_proteins_long = da_rows[, c("uniprot_acc", "comparison", "raw_pvalue", "fdr_qvalue", "log2FC")],
    list_of_volcano_plots = data.frame(title = "A_vs_B", stringsAsFactors = FALSE),
    list_of_volcano_plots_with_gene_names = tibble::tibble(title = character(), plot = list()),
    num_sig_da_molecules = data.frame(comparison = "A_vs_B", status = "Up", counts = 1),
    num_sig_da_genes_barplot_only_significant = empty_plot,
    num_of_comparison_only_significant = 1,
    num_sig_da_genes_barplot_with_not_significant = empty_plot,
    num_of_comparison_with_not_significant = 1
  )
  da_analysis_results$list_of_volcano_plots$plot <- list(empty_plot)

  withProtDaNamespaceBindings(
    list(
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        if (param_name_string %in% names(theObject@args)) {
          return(theObject@args[[param_name_string]])
        }
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      },
      ggsave = function(filename, plot, ...) {
        dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
        writeLines("plot", filename)
        captured$plot_names <<- c(captured$plot_names, basename(filename))
        invisible(filename)
      },
      plotSA = function(...) {
        plot.new()
      },
      savePlot = function(plot, base_path, plot_name, ...) {
        dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
        captured$plot_names <<- c(captured$plot_names, as.character(plot_name)[1])
        invisible(file.path(base_path, paste0(as.character(plot_name)[1], ".mock")))
      },
      createDirIfNotExists = function(path) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
      },
      writeInteractiveVolcanoPlotProteomics = function(..., counts_tbl, groups) {
        captured$interactive <<- list(counts_tbl = counts_tbl, groups = groups)
        TRUE
      }
    ),
    {
      testthat::local_mocked_bindings(
        vroom_write = function(x, file = NULL, path = NULL, ...) {
          output_path <- if (!is.null(path)) path else file
          dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
          utils::write.table(x, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
          captured$vroom_paths <<- c(captured$vroom_paths, basename(output_path))
          invisible(x)
        },
        .package = "vroom"
      )
      testthat::local_mocked_bindings(
        write_xlsx = function(x, path, ...) {
          dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
          writeLines("xlsx", path)
          captured$xlsx_paths <<- c(captured$xlsx_paths, basename(path))
          invisible(path)
        },
        .package = "writexl"
      )

      expect_true(output_da_analysis_results(
        da_analysis_results_list = da_analysis_results,
        theObject = obj,
        uniprot_tbl = obj@args$uniprot_tbl
      ))

      save_da_protein_list(
        list("A_vs_B" = data.frame(fdr_qvalue = c(0.2, 0.01), value = c("late", "early"))),
        row_id = "uniprot_acc",
        results_dir = output_dir,
        file_suffix = "_da.tsv"
      )
    }
  )

  expect_true("lfc_qval_long.tsv" %in% captured$vroom_paths)
  expect_true("da_proteins_wide.tsv" %in% captured$vroom_paths)
  expect_true("da_proteins_long_annot.tsv" %in% captured$vroom_paths)
  expect_true("A_vs_B_da.tsv" %in% captured$vroom_paths)
  expect_true("da_proteins_wide.xlsx" %in% captured$xlsx_paths)
  expect_equal(colnames(captured$interactive$counts_tbl), c("S1", "S2", "S3", "S4"))
  expect_equal(unname(captured$interactive$groups), c("A", "A", "B", "B"))
})

test_that("the DA base model takes its block from the declared technical replicate column and records why (KIT-035)", {
  fit_base <- getFromNamespace("fitProtDaBaseModel", "MultiScholaR")
  calls <- new.env(parent = emptyenv())
  mocks <- list(
    duplicateCorrelationFn = function(object, design, block, ...) { calls$dc <- (calls$dc %||% 0L) + 1L; list(consensus.correlation = 0.42) },
    lmFitFn = function(object, design, block = NULL, correlation = NULL, ...) list(block = block, correlation = correlation))
  data <- matrix(stats::rnorm(32), nrow = 4, dimnames = list(paste0("P", 1:4), paste0("S", 1:8)))
  design_m <- stats::model.matrix(~ 0 + group, data.frame(group = rep(c("A", "B"), each = 4)))
  rownames(design_m) <- colnames(data)
  base <- data.frame(Run = colnames(data), group = rep(c("A", "B"), each = 4), stringsAsFactors = FALSE, row.names = colnames(data))
  run <- function(design, technical_replicate_id = "tech_rep_group") do.call(fit_base, c(list(data, design_m, design, colnames(data)),
    mocks, list(group_id = "group", technical_replicate_id = technical_replicate_id)))
  # True technical replicates: two blocks per group, each with two runs; the declared values are already nested in the group.
  technical <- base; technical$tech_rep_group <- c("A_1", "A_1", "A_2", "A_2", "B_3", "B_3", "B_4", "B_4")
  fit <- run(technical)
  expect_true(fit$uses_duplicate_correlation)
  expect_identical(fit$block_reason, "technical_replicates")
  expect_identical(fit$block, technical$tech_rep_group)
  expect_identical(fit$technical_replicate_id, "tech_rep_group")
  expect_equal(fit$consensus_correlation, 0.42)
  expect_identical(fit$fit$correlation, 0.42)
  expect_identical(calls$dc, 1L)
  # A plain replicate index that recurs across groups is prefixed by the group, the block the model used before.
  index <- base; index$tech_rep_group <- c(1, 1, 2, 2, 1, 1, 2, 2)
  expect_identical(run(index)$block, c("A_1", "A_1", "A_2", "A_2", "B_1", "B_1", "B_2", "B_2"))
  expect_identical(calls$dc, 2L) # two blocks per group: technical replicates, duplicateCorrelation ran again
  # A block equal to the group is spanned by the design: standard lmFit, the reason recorded, no duplicateCorrelation.
  spanned <- base; spanned$tech_rep_group <- c("A_1", "A_1", "A_1", "A_1", "B_1", "B_1", "B_1", "B_1")
  fit <- run(spanned)
  expect_false(fit$uses_duplicate_correlation)
  expect_identical(fit$block_reason, "spanned_by_design")
  expect_null(fit$fit$block)
  expect_identical(calls$dc, 2L)
  # No duplicates, no declaration, missing values.
  unique <- base; unique$tech_rep_group <- paste0("r", 1:8)
  expect_identical(run(unique)$block_reason, "no_duplicates")
  expect_identical(run(base)$block_reason, "no_declaration")
  missing <- technical; missing$tech_rep_group[2] <- NA
  expect_identical(run(missing)$block_reason, "missing_values")
  expect_identical(calls$dc, 2L)
})

test_that("the DA helper fits limpa's dpcDE with the W covariates and the resolved block when the object carries the quantified EList (KIT-053)", {
  skip_if_not_installed("limpa")
  skip_if_not_installed("limma")
  ids <- sprintf("P%02d", 1:30)
  samples <- paste0("S", 1:8)
  set.seed(53)
  E <- matrix(rnorm(240, 20, 1), 30, 8, dimnames = list(ids, samples))
  E[1:5, 5:8] <- E[1:5, 5:8] + 2
  nobs <- matrix(2L, 30, 8, dimnames = dimnames(E)); nobs[sample(240, 20)] <- 0L
  elist <- methods::new("EList", list(E = E, genes = data.frame(NPrec = rep(2L, 30L), PropObs = rep(1, 30L), row.names = ids),
    other = list(standard.error = matrix(0.2, 30, 8, dimnames = dimnames(E)), n.observations = nobs)))
  design <- data.frame(Run = samples, group = rep(c("control", "treated"), each = 4L), replicates = rep(c(1L, 1L, 2L, 2L), 2L),
    W1 = seq(-0.5, 0.5, length.out = 8), stringsAsFactors = FALSE)
  table <- data.frame(Protein.Ids = ids, as.data.frame(E, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(table) <- NULL
  object <- ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = ids, stringsAsFactors = FALSE), design_matrix = design,
    sample_id = "Run", group_id = "group", technical_replicate_id = "replicates",
    args = list(limpa_dpc_quant_results = list(quantified_elist = elist, standard_errors = elist$other$standard.error)))
  captured <- new.env(parent = emptyenv())
  fake_dpcDE <- function(y, design, plot = TRUE, ...) {
    dots <- list(...)
    captured$block <- dots$block
    captured$design <- design
    captured$rows <- nrow(y)
    limma::lmFit(y$E, design)
  }
  testthat::local_mocked_bindings(dpcDE = fake_dpcDE, .package = "limpa")
  result <- suppressMessages(differentialAbundanceAnalysis(object,
    contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0.5, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
  # The branch ran: the formula gained the W term, the block reached dpcDE (the fixture declares technical replicates:
  # two runs share group and replicate number), TREAT was not applied, and the model is named.
  expect_identical(result$da_model$model, "dpcDE")
  expect_identical(result$da_model$formula_used, "~ 0 + group + W1")
  expect_identical(result$da_model$w_terms, "W1")
  expect_identical(result$da_model$k, 1L)
  expect_false(result$da_model$treat_applied)
  expect_identical(result$da_model$treat_lfc_cutoff_requested, 0.5)
  expect_identical(result$da_model$block_reason, "technical_replicates")
  expect_identical(captured$block, paste(design$group, design$replicates, sep = "_"))
  expect_identical(colnames(captured$design), c("groupcontrol", "grouptreated", "W1"))
  expect_identical(captured$rows, 30L)
  expect_identical(nrow(result$da_proteins_long), 30L)
  expect_true(all(c("Protein.Ids", "comparison", "raw_pvalue", "fdr_qvalue") %in% names(result$da_proteins_long)))
  expect_setequal(result$da_proteins_long$Protein.Ids, ids)
  expect_true(all(c("Protein.Ids", "comparison", "logFC", "raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment") %in% names(result$contrasts_results_table)))
  # Without the EList the same call takes the limma path and names it.
  plain <- object; plain@args$limpa_dpc_quant_results <- NULL
  limma_result <- suppressMessages(differentialAbundanceAnalysis(plain,
    contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
  expect_identical(limma_result$da_model$model, "limma")
  expect_identical(limma_result$da_model$formula_used, "~ 0 + group")
  expect_identical(limma_result$da_model$k, 0L)
  expect_identical(result$da_model$ruv_mode, "covariates")
  expect_identical(limma_result$da_model$ruv_mode, "none")
  # KIT-056: a corrected object (the mode declared, no W terms) fits the plain formula on the corrected EList and is
  # named; a corrected object whose design still carries W terms is refused (a double correction).
  corrected <- object
  corrected@design_matrix$W1 <- NULL
  corrected@args$ruvCovariates <- list(k = 1L, ruv_mode = "correction")
  corrected_result <- suppressMessages(differentialAbundanceAnalysis(corrected,
    contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0.5, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
  expect_identical(corrected_result$da_model$model, "dpcDE")
  expect_identical(corrected_result$da_model$formula_used, "~ 0 + group")
  expect_identical(corrected_result$da_model$ruv_mode, "correction")
  expect_identical(corrected_result$da_model$k, 0L)
  expect_identical(colnames(captured$design), c("groupcontrol", "grouptreated"))
  double <- object
  double@args$ruvCovariates <- list(k = 1L, ruv_mode = "correction")
  expect_error(suppressMessages(differentialAbundanceAnalysis(double,
    contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0.5, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")), "double correction")
})
