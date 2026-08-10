library(testthat)

makeQualityFunctionWithOverrides <- function(fun, replacements) {
  overridden_fun <- fun
  environment(overridden_fun) <- list2env(replacements, parent = environment(fun))
  overridden_fun
}

if (!methods::isClass("QualityProteinObject")) {
  methods::setClass(
    "QualityProteinObject",
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

test_that("getProteinsHeatMap calls ComplexHeatmap with supported arguments", {
  heatmap_fun <- makeQualityFunctionWithOverrides(
    getProteinsHeatMap,
    list(Legend = ComplexHeatmap::Legend)
  )
  protein_matrix <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S2", "S1"))
  )
  metadata_tbl <- data.frame(
    Run = c("S1", "S2", "S3"),
    is_HEK = c(FALSE, FALSE, TRUE),
    condition = c("A", "B", "A"),
    batch = c("X", "Y", "X"),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(suppressWarnings(heatmap_fun(
    protein_matrix = protein_matrix,
    metadata_tbl = metadata_tbl,
    is_HEK_column = is_HEK,
    metadata_column_selected = c("condition", "batch"),
    metadata_column_labels = c("Condition", "Batch"),
    colour_rules = list(
      Condition = c(A = "#1b9e77", B = "#d95f02"),
      Batch = c(X = "#7570b3", Y = "#e7298a")
    ),
    columns_to_exclude = character(),
    core_utilisation_samples = FALSE,
    sample_id_column = Run,
    raster_device = "CairoPNG"
  )))

  expect_s4_class(result$heatmap, "Heatmap")
  expect_length(result$legend, 2L)
})

test_that("outputDeAnalysisResults forwards the DA q-value threshold", {
  captured <- new.env(parent = emptyenv())
  captured$threshold <- NULL

  the_object <- methods::new(
    "QualityProteinObject",
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
      stringsAsFactors = FALSE
    )
  )

  de_output_dir <- tempfile("quality-da-output-")
  publication_graphs_dir <- tempfile("quality-da-pubs-")
  dir.create(de_output_dir, recursive = TRUE)
  dir.create(publication_graphs_dir, recursive = TRUE)
  withr::defer({
    unlink(de_output_dir, recursive = TRUE, force = TRUE)
    unlink(publication_graphs_dir, recursive = TRUE, force = TRUE)
  })

  de_proteins_long <- data.frame(
    uniprot_acc = c("P1", "P2"),
    comparison = c("A_vs_B", "A_vs_B"),
    fdr_qvalue = c(0.01, 0.02),
    raw_pvalue = c(0.001, 0.01),
    log2FC = c(2, -1.5),
    stringsAsFactors = FALSE
  )
  analysis_results <- list(
    theObject = the_object,
    contrasts_results = list(fit.eb = list(coefficients = matrix(1, ncol = 1))),
    de_proteins_long = de_proteins_long
  )
  uniprot_tbl <- data.frame(
    Entry = c("P1", "P2"),
    stringsAsFactors = FALSE
  )

  output_fun <- makeQualityFunctionWithOverrides(
    outputDeAnalysisResults,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, key, default = NULL) {
        switch(key,
          uniprot_tbl = uniprot_tbl,
          de_output_dir = de_output_dir,
          publication_graphs_dir = publication_graphs_dir,
          file_prefix = "de_proteins",
          plots_format = character(),
          args_row_id = "uniprot_acc",
          de_q_val_thresh = 0.0123,
          gene_names_column = "gene_name",
          fdr_column = "fdr_qvalue",
          raw_p_value_column = "raw_pvalue",
          log2fc_column = "log2FC",
          uniprot_id_column = "Entry",
          display_columns = "uniprot_acc",
          default
        )
      },
      updateParamInObject = function(theObject, key) theObject,
      writeInteractiveVolcanoPlotProteomics = function(
        da_proteins_long,
        uniprot_tbl,
        fit.eb,
        publication_graphs_dir,
        args_row_id,
        fdr_column,
        raw_p_value_column,
        log2fc_column,
        da_q_val_thresh,
        counts_tbl,
        groups,
        uniprot_id_column,
        gene_names_column,
        display_columns
      ) {
        captured$threshold <- da_q_val_thresh
        invisible(NULL)
      }
    )
  )

  output_fun(
    de_analysis_results_list = analysis_results,
    theObject = the_object,
    uniprot_tbl = uniprot_tbl,
    de_output_dir = de_output_dir,
    publication_graphs_dir = publication_graphs_dir
  )

  expect_identical(captured$threshold, 0.0123)
})
