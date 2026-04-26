# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

if (!methods::isClass("MockProtCompatObject")) {
  methods::setClass(
    "MockProtCompatObject",
    slots = c(
      protein_quant_table = "data.frame",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      protein_id_column = "character",
      protein_id_table = "data.frame"
    )
  )
}

buildCompatPlot <- function(title = "plot") {
  ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle(title)
}

test_that("protein_deAnalysisWrapperFunction preserves orchestration and DA result assembly", {
  wrapper_fun <- makeFunctionWithOverrides(
    protein_deAnalysisWrapperFunction,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, key, default = NULL) {
        switch(key,
          contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
          formula_string = "~ 0 + group",
          group_id = "group",
          de_q_val_thresh = 0.05,
          treat_lfc_cutoff = 0,
          eBayes_trend = TRUE,
          eBayes_robust = TRUE,
          args_group_pattern = "(\\d+)",
          args_row_id = "uniprot_acc",
          default
        )
      },
      updateParamInObject = function(theObject, key) theObject,
      plotRle = function(theObject, group_id) buildCompatPlot("rle"),
      plotPca = function(theObject, grouping_variable, label_column, title, font_size) buildCompatPlot("pca"),
      plotNumOfValuesNoLog = function(tbl) buildCompatPlot("values"),
      pull = dplyr::pull,
      one_of = tidyselect::one_of,
      column_to_rownames = tibble::column_to_rownames,
      set_colnames = function(x, nm) {
        colnames(x) <- nm
        x
      },
      rownames_to_column = tibble::rownames_to_column,
      runTestsContrasts = function(mat,
                                   contrast_strings,
                                   design_matrix,
                                   formula_string,
                                   weights,
                                   treat_lfc_cutoff,
                                   eBayes_trend,
                                   eBayes_robust) {
        list(
          results = data.frame(
            uniprot_acc = c("P1", "P2"),
            logFC = c(2, -1.5),
            raw_pvalue = c(0.001, 0.01),
            fdr_qvalue = c(0.01, 0.02),
            comparison = c("A_vs_B", "A_vs_B"),
            log_intensity = c(10, 12),
            analysis_type = c("RUV applied", "RUV applied"),
            stringsAsFactors = FALSE
          ),
          fit.eb = list(coefficients = matrix(1, ncol = 1))
        )
      },
      getSignificantData = function(...) {
        data.frame(
          uniprot_acc = c("P1", "P2"),
          raw_pvalue = c(0.001, 0.01),
          fdr_qvalue = c(0.01, 0.02),
          lqm = c(2, 1.7),
          logFC = c(2, -1.5),
          comparison = c("A_vs_B", "A_vs_B"),
          log_intensity = c(10, 12),
          analysis_type = c("RUV applied", "RUV applied"),
          colour = c("purple", "purple"),
          stringsAsFactors = FALSE
        )
      },
      plotVolcano = function(...) buildCompatPlot("volcano"),
      printCountDaGenesTable = function(...) {
        list(table = data.frame(comparison = "A_vs_B", status = "Significant", counts = 2L, stringsAsFactors = FALSE))
      },
      printPValuesDistribution = function(...) buildCompatPlot("pvalhist"),
      createDaResultsLongFormat = function(...) {
        data.frame(
          uniprot_acc = c("P1", "P2"),
          comparison = c("A_vs_B", "A_vs_B"),
          fdr_qvalue = c(0.01, 0.02),
          raw_pvalue = c(0.001, 0.01),
          log2FC = c(2, -1.5),
          analysis_type = c("RUV applied", "RUV applied"),
          stringsAsFactors = FALSE
        )
      },
      plotOneVolcanoNoVerticalLines = function(x, y, ...) buildCompatPlot(y),
      raw_pvalue = rlang::sym("raw_pvalue"),
      fdr_qvalue = rlang::sym("fdr_qvalue"),
      lqm = rlang::sym("lqm"),
      logFC = rlang::sym("logFC"),
      analysis_type = rlang::sym("analysis_type"),
      log_intensity = rlang::sym("log_intensity"),
      sym = rlang::sym
    )
  )

  the_object <- methods::new(
    "MockProtCompatObject",
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      S1 = c(10, 11),
      S2 = c(12, 13),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      sample = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "sample",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      gene_name = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    )
  )
  rownames(the_object@protein_quant_table) <- the_object@protein_quant_table$uniprot_acc

  results <- wrapper_fun(the_object)

  expect_true(all(c(
    "theObject",
    "rle_plot",
    "pca_plot",
    "pca_plot_with_labels",
    "plot_num_of_values",
    "contrasts_results",
    "contrasts_results_table",
    "significant_rows",
    "volplot_plot",
    "num_sig_de_molecules_first_go",
    "pvalhist",
    "norm_counts",
    "de_proteins_wide",
    "de_proteins_long",
    "list_of_volcano_plots",
    "num_sig_de_molecules"
  ) %in% names(results)))
  expect_s3_class(results$rle_plot, "ggplot")
  expect_s3_class(results$pca_plot, "ggplot")
  expect_equal(nrow(results$de_proteins_long), 2L)
  expect_equal(results$norm_counts$uniprot_acc, c("P1", "P2"))
  expect_identical(results$list_of_volcano_plots$title[[1]], "A_vs_B")
})

test_that("outputDeAnalysisResults preserves artifact routing and interactive writer delegation", {
  captured <- new.env(parent = emptyenv())
  captured$ggsave <- character()
  captured$interactive <- NULL

  the_object <- methods::new(
    "MockProtCompatObject",
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      S1 = c(10, 11),
      S2 = c(12, 13),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      sample = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "sample",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      gene_name = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    )
  )

  de_output_dir <- tempfile("prot-da-output-")
  publication_graphs_dir <- tempfile("prot-da-pubs-")
  dir.create(de_output_dir, recursive = TRUE)
  dir.create(publication_graphs_dir, recursive = TRUE)
  withr::defer({
    unlink(de_output_dir, recursive = TRUE, force = TRUE)
    unlink(publication_graphs_dir, recursive = TRUE, force = TRUE)
  })

  de_analysis_results_list <- list(
    pca_plot = buildCompatPlot("pca"),
    pca_plot_with_labels = buildCompatPlot("pca-labeled"),
    rle_plot = buildCompatPlot("rle"),
    plot_num_of_values = buildCompatPlot("values"),
    contrasts_results = list(fit.eb = list(coefficients = matrix(1, ncol = 1))),
    significant_rows = data.frame(
      uniprot_acc = c("P1", "P2"),
      raw_pvalue = c(0.001, 0.01),
      fdr_qvalue = c(0.01, 0.02),
      log2FC = c(2, -1.5),
      comparison = c("A_vs_B", "A_vs_B"),
      stringsAsFactors = FALSE
    ),
    volplot_plot = buildCompatPlot("volcano"),
    num_sig_de_molecules_first_go = list(
      table = data.frame(comparison = "A_vs_B", counts = 2L, stringsAsFactors = FALSE)
    ),
    num_sig_de_molecules = data.frame(
      comparison = "A_vs_B",
      status = "Significant and Up",
      counts = 2L,
      stringsAsFactors = FALSE
    ),
    num_sig_de_genes_barplot_only_significant = buildCompatPlot("num-sig"),
    num_of_comparison_only_significant = 1L,
    pvalhist = buildCompatPlot("pvalhist"),
    de_proteins_wide = data.frame(
      uniprot_acc = c("P1", "P2"),
      `log2FC:A_vs_B` = c(2, -1.5),
      stringsAsFactors = FALSE
    ),
    de_proteins_long = data.frame(
      uniprot_acc = c("P1", "P2"),
      comparison = c("A_vs_B", "A_vs_B"),
      fdr_qvalue = c(0.01, 0.02),
      raw_pvalue = c(0.001, 0.01),
      log2FC = c(2, -1.5),
      stringsAsFactors = FALSE
    ),
    list_of_volcano_plots = tibble::tibble(
      title = "A_vs_B",
      plot = list(buildCompatPlot("static-volcano"))
    ),
    theObject = the_object
  )
  uniprot_tbl <- data.frame(
    Entry = c("P1", "P2"),
    best_uniprot_acc = c("P1", "P2"),
    gene_name = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )

  output_fun <- makeFunctionWithOverrides(
    outputDeAnalysisResults,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, key, default = NULL) {
        switch(key,
          uniprot_tbl = uniprot_tbl,
          de_output_dir = de_output_dir,
          publication_graphs_dir = publication_graphs_dir,
          file_prefix = "de_proteins",
          plots_format = c("pdf", "png"),
          args_row_id = "uniprot_acc",
          de_q_val_thresh = 0.05,
          gene_names_column = "gene_name",
          fdr_column = "fdr_qvalue",
          raw_p_value_column = "raw_pvalue",
          log2fc_column = "log2FC",
          uniprot_id_column = "Entry",
          display_columns = "best_uniprot_acc",
          default
        )
      },
      updateParamInObject = function(theObject, key) theObject,
      plotSA = function(fit) invisible(fit),
      ggsave = function(filename, plot, ...) {
        captured$ggsave <- c(captured$ggsave, filename)
        invisible(filename)
      },
      pdf = function(...) invisible(NULL),
      png = function(...) invisible(NULL),
      dev.off = function(...) invisible(NULL),
      str_split = stringr::str_split,
      sym = rlang::sym,
      column_to_rownames = tibble::column_to_rownames,
      createDirIfNotExists = function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE),
      writeInteractiveVolcanoPlotProteomics = function(da_proteins_long,
                                                       uniprot_tbl,
                                                       publication_graphs_dir,
                                                       ...) {
        captured$interactive <- list(
          rows = nrow(da_proteins_long),
          publication_graphs_dir = publication_graphs_dir,
          args = list(...)
        )
        invisible(NULL)
      }
    )
  )

  output_fun(
    de_analysis_results_list = de_analysis_results_list,
    theObject = the_object,
    uniprot_tbl = uniprot_tbl,
    de_output_dir = de_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    file_prefix = "de_proteins",
    plots_format = c("pdf", "png"),
    args_row_id = "uniprot_acc",
    display_columns = "best_uniprot_acc"
  )

  expect_true(any(grepl("PCA_plot.pdf", captured$ggsave, fixed = TRUE)))
  expect_true(any(grepl("volplot_gg_all.png", captured$ggsave, fixed = TRUE)))
  expect_true(any(grepl("num_sig_de_genes_barplot.pdf", captured$ggsave, fixed = TRUE)))
  expect_true(file.exists(file.path(de_output_dir, "lfc_qval_long.tsv")))
  expect_true(file.exists(file.path(de_output_dir, "de_proteins_wide.tsv")))
  expect_true(file.exists(file.path(de_output_dir, "de_proteins_long_annot.xlsx")))
  expect_true(file.exists(file.path(de_output_dir, "fit.eb.RDS")))
  expect_true(file.exists(file.path(publication_graphs_dir, "NumSigDeMolecules", "num_sig_de_molecules.tab")))
  expect_identical(captured$interactive$rows, 2L)
  expect_identical(captured$interactive$publication_graphs_dir, publication_graphs_dir)
})
